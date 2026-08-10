"""duckdb -> per-species feature matrices for the BFD comparative dashboard.

All matrices are indexed by LOCUSTAG (one row per genome in the current subset).
No I/O happens at import time; call get_connection() to open the db.

This module is deliberately dataset-agnostic: every function takes the db path /
directory it needs as an argument (or a ``rows`` LOCUSTAG subset), so a build is
driven entirely by a DashboardConfig (lib/config.py) instead of module-level
constants that encode one project's layout.

Zero-handling rule (documented once, used everywhere a family/species cell can be
a true zero -- Pfam/CAZy/MEROPS domain counts): a family absent from a species is
a true zero, not missing data. Where a method requires strictly positive values
(chi-square/CA distances, log/fold-change), add PSEUDOCOUNT once at matrix-build
time rather than letting each downstream consumer improvise its own fix.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd

PSEUDOCOUNT = 0.5  # Anscombe-style; applied once, at matrix build, not per-consumer

# phyling names this file "protein-{dataset}-taxa_{n}.{lineage}.fasttree.support.treefile"
# -- both {dataset} and {n} change every time the genome set is edited or renamed, so
# find_tree_path() globs for it instead of pinning a literal path.
TREE_GLOB = "phyling_pep/protein/buildtree/*/fasttree/*.fasttree.support.treefile"

# Fallback used only when a config does not set merops_families_path and this file
# happens to exist on the host. The MEROPS install path is lab/version specific, so
# treat it as a convenience, never a requirement.
_DEFAULT_MEROPS_FAMILIES_PATH = Path("/srv/projects/db/MEROPS/124/merops_lib.families.tab")


# --------------------------------------------------------------------------- connection / schema
def get_connection(db_path: Path | str) -> duckdb.DuckDBPyConnection:
    return duckdb.connect(str(db_path), read_only=True)


def list_tables(con: duckdb.DuckDBPyConnection) -> list[str]:
    return [r[0] for r in con.execute("SHOW TABLES").fetchall()]


def describe_table(con: duckdb.DuckDBPyConnection, table: str) -> pd.DataFrame:
    return con.execute(f"DESCRIBE {table}").fetchdf()


def feature_coverage(con: duckdb.DuckDBPyConnection) -> dict[str, int]:
    """How many species have data in each analysis table.

    A large BFD database may have assembly/BUSCO data for every genome while the
    functional-annotation tables (pfam/cazy/merops/codon/...) are only populated
    for a subset that has finished the pipeline. This lets the builder tell the
    user exactly which views will be empty for a chosen subset.
    """
    tables = ["asm_stats", "busco_genome", "codon_frequency", "aa_frequency",
              "pfam", "cazy", "merops", "signalp", "tmhmm", "targetp", "idp",
              "wolfpsort", "predgpi"]
    # asm/busco join on LOCUSTAG; the functional tables all use species_prefix
    # (which equals LOCUSTAG for the same genome -- see functional_species).
    locustag_keyed = {"asm_stats", "busco_genome"}
    out = {}
    for t in tables:
        if t not in list_tables(con):
            continue
        keycol = "LOCUSTAG" if t in locustag_keyed else "species_prefix"
        try:
            out[t] = int(
                con.execute(
                    f"SELECT COUNT(DISTINCT {keycol}) FROM {t}"
                ).fetchone()[0]
            )
        except duckdb.Error:
            out[t] = 0
    return out


# --------------------------------------------------------------------------- taxa filtering
def resolve_taxa_filter(
    con: duckdb.DuckDBPyConnection,
    genus: str | None = None,
    species: str | None = None,
    locustags: list[str] | None = None,
    exclude_locustags: list[str] | None = None,
) -> list[str]:
    """Resolve a taxon description into a concrete list of LOCUSTAGs, for consistent
    subsetting across every FEATURE_SETS loader and analysis tool. `genus`/`species`
    do a case-insensitive substring match against the `species` table; `locustags`
    (if given) intersects with that match instead of replacing it, so callers can
    combine e.g. genus="Exophiala" with an explicit locustags allowlist. Returns
    LOCUSTAGs in species-table order (stable across calls).
    """
    sp = con.execute("SELECT LOCUSTAG, GENUS, SPECIES FROM species ORDER BY LOCUSTAG").fetchdf()
    mask = pd.Series(True, index=sp.index)
    if genus:
        mask &= sp["GENUS"].str.contains(genus, case=False, na=False)
    if species:
        mask &= sp["SPECIES"].str.contains(species, case=False, na=False)
    result = sp.loc[mask, "LOCUSTAG"].tolist()
    if locustags:
        wanted = set(locustags)
        result = [lt for lt in result if lt in wanted]
    if exclude_locustags:
        excluded = set(exclude_locustags)
        result = [lt for lt in result if lt not in excluded]
    return result


def functional_species(con: duckdb.DuckDBPyConnection) -> set[str]:
    """LOCUSTAGs present in at least one functional-annotation hit table. These are
    the genomes whose proteomes have actually been screened -- a large BFD db may
    only have this for a fraction of `species` rows."""
    hit_tables = ["pfam", "cazy", "merops", "signalp", "tmhmm", "targetp", "idp",
                  "wolfpsort", "predgpi", "codon_frequency", "aa_frequency"]
    found = set()
    for t in hit_tables:
        if t in list_tables(con):
            col = "species_prefix" if t in {"codon_frequency", "aa_frequency"} else "species_prefix"
            try:
                found.update(
                    r[0] for r in con.execute(
                        f"SELECT DISTINCT species_prefix FROM {t}"
                    ).fetchall()
                )
            except duckdb.Error:
                pass
    return found


# --------------------------------------------------------------------------- species / overview
def get_species_table(
    con: duckdb.DuckDBPyConnection, rows: list[str] | None = None
) -> pd.DataFrame:
    """Per-genome overview table, indexed by LOCUSTAG.

    Backbone is the `species` taxonomy table (every genome in the database), enriched
    with assembly stats from `asm_stats` (present for nearly all genomes) and the
    richer per-genome annotation summary `v_species_summary` (gene count, BUSCO-adjacent
    counts, etc. -- only populated for genomes that have finished the pipeline). Stats
    are pulled from asm_stats whenever the summary is missing, so a large BFD database
    whose summary view is still sparse produces real genome stats for every taxon
    instead of dropping genomes.

    Overlapping columns (TOTAL_LENGTH, GC_PERCENT, N50_bp, masked_pct) prefer the
    v_species_summary value when present, else the asm_stats value."""
    tables = set(list_tables(con))

    sp = con.execute(
        'SELECT LOCUSTAG, ASMID, PHYLUM, SUBPHYLUM, CLASS, "ORDER", FAMILY, GENUS, '
        "SPECIES, STRAIN, NCBI_TAXONID, BUSCO_LINEAGE FROM species ORDER BY LOCUSTAG"
    ).fetchdf()

    summary = None
    if "v_species_summary" in tables:
        summary = con.execute("SELECT * FROM v_species_summary").fetchdf().set_index("LOCUSTAG")

    asm = None
    if "asm_stats" in tables:
        asm = con.execute(
            "SELECT * FROM asm_stats ORDER BY LOCUSTAG"
        ).fetchdf().set_index("LOCUSTAG")
        asm = asm.drop(columns=["ASMID"], errors="ignore")

    df = sp.set_index("LOCUSTAG")
    if asm is not None:
        # asm_stats brings real genome stats for every genome, even those whose
        # v_species_summary row has not been built yet. Drop any column that would
        # collide with the species backbone (e.g. ASMID).
        asm = asm.drop(columns=[c for c in asm.columns if c in df.columns], errors="ignore")
        df = df.join(asm, how="left")
    if summary is not None:
        # v_species_summary re-states the taxonomy (SPECIES/GENUS/...) -- the backbone
        # already has that for all genomes, so only merge the summary-only annotation
        # columns (gene_count, mean_protein_length, pfam_proteins, ...).
        summary_keep = [c for c in summary.columns if c not in df.columns]
        df = df.join(summary[summary_keep], how="left")
        # For the columns asm_stats and v_species_summary share (TOTAL_LENGTH,
        # GC_PERCENT, N50_bp, masked_pct, contig_count), prefer the summary value
        # when the genome has one, else keep the asm_stats value. combine_first
        # aligns on the full df index, so the 1-row summary can't clobber the rest.
        if asm is not None:
            shared = [c for c in summary.columns if c in df.columns and c in asm.columns]
            for c in shared:
                df[c] = summary[c].reindex(df.index).combine_first(df[c])
    df = df.reset_index().sort_values("LOCUSTAG").set_index("LOCUSTAG")
    return _filter_rows(df, rows)


def get_busco_scores(busco_dir: Path | None) -> pd.DataFrame:
    """Parse per-genome protein-set BUSCO summaries into a LOCUSTAG-indexed table.

    Returns NaN-filled rows for any species missing a summary file rather than
    dropping them -- callers should show "not available" in the UI. A None/absent
    directory yields an empty frame (all species "not available")."""
    rows = {}
    if busco_dir and busco_dir.exists():
        for txt in busco_dir.glob("*/*.BUSCO_summary.*.txt"):
            locustag = txt.name.split(".BUSCO_summary.")[0]
            m = _BUSCO_HEADER_RE.search(txt.read_text())
            if not m:
                continue
            c, s, d, f, miss, n = m.groups()
            rows[locustag] = {
                "busco_complete_pct": float(c),
                "busco_single_pct": float(s),
                "busco_duplicated_pct": float(d),
                "busco_fragmented_pct": float(f),
                "busco_missing_pct": float(miss),
                "busco_n": int(n),
            }
    return pd.DataFrame.from_dict(rows, orient="index")


_BUSCO_HEADER_RE = re.compile(
    r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%,n:(\d+)"
)


def _filter_rows(df: pd.DataFrame, rows: list[str] | None) -> pd.DataFrame:
    """Subset a genome-indexed matrix to `rows` (LOCUSTAGs/species_prefixes) BEFORE
    any normalization that follows -- z-scoring/CA computed on the full cohort and
    then sliced would leave a subset's "standardized" values centered/scaled on
    genomes not even in the subset, which is invalid. Every FEATURE_SETS loader
    takes `rows` for exactly this reason."""
    if rows is None:
        return df
    return df.loc[df.index.isin(rows)]


# --------------------------------------------------------------------------- feature matrices
def get_assembly_matrix(con: duckdb.DuckDBPyConnection, rows: list[str] | None = None) -> pd.DataFrame:
    """Numeric assembly-stats block, standardized (z-score) per column. `rows`
    (LOCUSTAGs) subsets BEFORE z-scoring, so a taxon-filtered ordination is
    standardized on that subset, not silently inheriting the full-cohort mean/std."""
    df = con.execute("SELECT * FROM asm_stats ORDER BY LOCUSTAG").fetchdf().set_index("LOCUSTAG")
    df = df.drop(columns=["ASMID"], errors="ignore")
    df = _filter_rows(df, rows)
    numeric = df.select_dtypes(include="number")
    z = (numeric - numeric.mean()) / numeric.std(ddof=0).replace(0, 1)
    return z


def get_codon_matrix(con: duckdb.DuckDBPyConnection, rows: list[str] | None = None) -> pd.DataFrame:
    """species x codon frequency profile (rows ~sum to 1). Raw frequencies -- CA in
    ordination.py handles the compositional/chi-square structure directly, no
    transform here."""
    df = con.execute("SELECT species_prefix, codon, frequency FROM codon_frequency").fetchdf()
    matrix = df.pivot(index="species_prefix", columns="codon", values="frequency").fillna(0.0)
    return _filter_rows(matrix, rows)


def get_aa_matrix(con: duckdb.DuckDBPyConnection, rows: list[str] | None = None) -> pd.DataFrame:
    df = con.execute("SELECT species_prefix, amino_acid, frequency FROM aa_frequency").fetchdf()
    matrix = df.pivot(index="species_prefix", columns="amino_acid", values="frequency").fillna(0.0)
    return _filter_rows(matrix, rows)


_DOMAIN_TABLES = {
    "pfam": ("pfam", "pfam_id"),
    "cazy": ("cazy", "HMM_id"),
    "merops": ("v_merops_family", "family"),
}


def get_pfam_accessions(con: duckdb.DuckDBPyConnection) -> dict[str, str]:
    """pfam_id (human-readable family name, e.g. "1-cysPrx_C") -> unversioned
    pfam_acc (e.g. "PF10417"), for linking out to the Pfam/InterPro entry page. A
    pfam_id can appear with more than one accession version across proteins (rare
    HMM updates); any one instance is fine since the URL only needs the PFxxxxx
    prefix, not the version suffix."""
    df = con.execute("SELECT DISTINCT pfam_id, pfam_acc FROM pfam").fetchdf()
    accession = df["pfam_acc"].str.split(".").str[0]
    return dict(zip(df["pfam_id"], accession))


class MeropsUnavailableError(RuntimeError):
    """Raised when a merops feature set is requested but the family-lookup file
    that builds v_merops_family is unavailable. Callers skip/disable the merops
    feature set instead of letting the whole build fail."""


def _ensure_merops_family_view(
    con: duckdb.DuckDBPyConnection, path: Path | None
) -> None:
    """(Re)create v_merops_family: the raw `merops` hits table joined to its MEROPS
    family code. merops_id there is the nearest-match reference sequence (e.g.
    MER0034723), which is ~1:1 with individual proteins -- grouping domain hits by
    merops_id, the way pfam groups by pfam_id, tells you almost nothing at the
    family level. merops_lib.families.tab (shipped with the local MEROPS database
    build) maps each merops_id to a family.member code such as "S09.996"; the
    family itself is the part before the dot ("S09"). Hits with no lookup match keep
    their raw merops_id as a fallback rather than being dropped.

    `path` is optional; if a path is not configured, a default well-known install
    is tried, and MeropsUnavailableError is raised if neither resolves."""
    if path is None:
        path = _DEFAULT_MEROPS_FAMILIES_PATH
    if not path.exists():
        raise MeropsUnavailableError(
            f"MEROPS family lookup file not found (set merops_families_path): {path}"
        )
    fam = pd.read_csv(path, sep="\t", header=None, names=["merops_id", "subfamily", "family_member"])
    fam["family"] = fam["family_member"].str.split(".").str[0]
    con.register("_merops_family_lookup", fam[["merops_id", "family"]])
    # TEMP (session-local, in the attached in-memory catalog) -- the db is opened
    # read_only, which forbids CREATE VIEW against its own catalog.
    con.execute(
        """
        CREATE OR REPLACE TEMP VIEW v_merops_family AS
        SELECT m.*, COALESCE(f.family, m.merops_id) AS family
        FROM merops m
        LEFT JOIN _merops_family_lookup f USING (merops_id)
        """
    )


def _domain_table(con: duckdb.DuckDBPyConnection, kind: str, merops_families_path: Path | None = None) -> tuple[str, str]:
    if kind == "merops":
        _ensure_merops_family_view(con, merops_families_path)
    return _DOMAIN_TABLES[kind]


def get_domain_matrix(
    con: duckdb.DuckDBPyConnection, kind: str, rows: list[str] | None = None,
    merops_families_path: Path | None = None,
) -> pd.DataFrame:
    """species x family RAW counts, gene-count-normalized (per 1000 genes) to avoid
    genome-size artifacts. Per-1000-genes normalization is row-wise (each genome
    divided by its own gene_count), so `rows` subsetting order doesn't change
    values here -- kept for API consistency with the other loaders."""
    table, family_col = _domain_table(con, kind, merops_families_path)
    counts = con.execute(
        f"SELECT species_prefix, {family_col} AS family, COUNT(*) AS n "
        f"FROM {table} GROUP BY species_prefix, {family_col}"
    ).fetchdf()
    matrix = counts.pivot(index="species_prefix", columns="family", values="n").fillna(0.0)
    matrix = _filter_rows(matrix, rows)

    gene_counts = get_species_table(con)["gene_count"]
    gene_counts = gene_counts.reindex(matrix.index)
    per_1000_genes = matrix.div(gene_counts, axis=0) * 1000.0
    return per_1000_genes


def get_domain_raw_counts(
    con: duckdb.DuckDBPyConnection, kind: str,
    merops_families_path: Path | None = None,
) -> pd.DataFrame:
    """Un-normalized companion to get_domain_matrix, for sanity-checking against
    gene_count."""
    table, family_col = _domain_table(con, kind, merops_families_path)
    counts = con.execute(
        f"SELECT species_prefix, {family_col} AS family, COUNT(*) AS n "
        f"FROM {table} GROUP BY species_prefix, {family_col}"
    ).fetchdf()
    return counts.pivot(index="species_prefix", columns="family", values="n").fillna(0.0)


def get_domain_example_proteins(
    con: duckdb.DuckDBPyConnection, kind: str, species_prefix: str, family: str,
    limit: int = 3, merops_families_path: Path | None = None,
) -> list[str]:
    table, family_col = _domain_table(con, kind, merops_families_path)
    df = con.execute(
        f"SELECT DISTINCT protein_id FROM {table} "
        f"WHERE species_prefix = ? AND {family_col} = ? LIMIT {int(limit)}",
        [species_prefix, family],
    ).fetchdf()
    return df["protein_id"].tolist()


MIN_REST_GROUP_SIZE = 3  # below this, "differs from the rest" is really "differs
# from 1-2 arbitrary other genomes," not a clade-level statement.


def differential_families(
    con: duckdb.DuckDBPyConnection, kind: str, top_n: int = 15,
    rows: list[str] | None = None, merops_families_path: Path | None = None,
) -> dict[str, dict]:
    """Per-species (focus) differential domain-family ranking against the
    gene-count-normalized matrix's mean-of-rest. NOT a z-score (a small 'rest'
    group is too small/non-normal for that to mean what it implies) -- ranks by
    (1) presence/absence flags as the headline list, (2) fold-change over group
    mean (with PSEUDOCOUNT) as a secondary continuous rank. Each entry carries a
    few example protein_ids for drill-down to specific genes.

    `rows` restricts which genomes participate (both as focus species AND as the
    "rest" pool). Every result reports `rest_group` (the actual LOCUSTAGs averaged)
    and `rest_group_n`, since a taxon-filtered "rest" can shrink below the
    MIN_REST_GROUP_SIZE floor and stop meaning "the clade" -- raises ValueError in
    that case rather than silently returning an over-broad claim."""
    normalized = get_domain_matrix(con, kind, rows=rows, merops_families_path=merops_families_path)
    if len(normalized) - 1 < MIN_REST_GROUP_SIZE:
        raise ValueError(
            f"Only {len(normalized)} genome(s) in this subset -- the 'rest of clade' "
            f"comparison group would have {len(normalized) - 1}, below the minimum of "
            f"{MIN_REST_GROUP_SIZE}. Widen the taxon filter before running differential_families."
        )
    pfam_accessions = get_pfam_accessions(con) if kind == "pfam" else {}
    result = {}
    for species in normalized.index:
        rest_group = [s for s in normalized.index if s != species]
        rest_mean = normalized.loc[rest_group].mean()
        mine = normalized.loc[species]

        present_only = mine[(mine > 0) & (rest_mean == 0)].sort_values(ascending=False)
        absent_only = rest_mean[(mine == 0) & (rest_mean > 0)].sort_values(ascending=False)
        fold = (mine + PSEUDOCOUNT) / (rest_mean + PSEUDOCOUNT)
        log2fold = np.log2(fold).sort_values(key=np.abs, ascending=False)

        def family_rows(series, value_label):
            out = []
            for family, value in series.head(top_n).items():
                row = {
                    "family": family,
                    value_label: float(value),
                    "species_value_per_1000": float(mine[family]),
                    "rest_mean_per_1000": float(rest_mean[family]),
                    "example_proteins": get_domain_example_proteins(
                        con, kind, species, family, merops_families_path=merops_families_path
                    ),
                }
                if kind == "pfam" and family in pfam_accessions:
                    row["pfam_accession"] = pfam_accessions[family]
                out.append(row)
            return out

        result[species] = {
            "rest_group": rest_group,
            "rest_group_n": len(rest_group),
            "present_only": family_rows(present_only, "value_per_1000"),
            "absent_only": family_rows(absent_only, "value_per_1000"),
            "top_fold_change": family_rows(log2fold, "log2_fold_change"),
        }
    return result


_LOCALIZATION_QUERIES = {
    "signalp_frac": (
        "SELECT species_prefix, COUNT(DISTINCT protein_id) AS n FROM signalp GROUP BY species_prefix",
        "n",
    ),
    "tmhmm_frac": (
        "SELECT species_prefix, COUNT(*) AS n FROM tmhmm WHERE PredHel > 0 GROUP BY species_prefix",
        "n",
    ),
    "secreted_frac": (
        "SELECT species_prefix, COUNT(*) AS n FROM targetp WHERE prediction = 'SP' GROUP BY species_prefix",
        "n",
    ),
    "gpi_frac": (
        "SELECT species_prefix, COUNT(DISTINCT protein_id) AS n FROM predgpi "
        "WHERE feature = 'GPI-anchor' GROUP BY species_prefix",
        "n",
    ),
    "has_disordered_region_frac": (
        "SELECT species_prefix, COUNT(DISTINCT protein_id) AS n FROM idp GROUP BY species_prefix",
        "n",
    ),
    "wolfpsort_extr_frac": (
        "SELECT species_prefix, COUNT(*) AS n FROM wolfpsort WHERE localization = 'extr' GROUP BY species_prefix",
        "n",
    ),
}


def get_localization_counts(con: duckdb.DuckDBPyConnection) -> pd.DataFrame:
    """Per-species raw integer counts for signal peptide (SignalP) / TM helix /
    secretion (TargetP) / GPI-anchor / high-disorder (IDP) proteins -- backs both
    get_localization_fractions (gene_count-normalized) and the composition
    payload's "counts" field, so bars can show "N proteins (X%)" instead of a
    bare fraction."""
    idx = get_species_table(con).index
    cols = {}
    for name, (sql, col) in _LOCALIZATION_QUERIES.items():
        series = con.execute(sql).fetchdf().set_index("species_prefix")[col]
        cols[name] = series.reindex(idx).fillna(0).astype(int)
    return pd.DataFrame(cols)


def get_localization_fractions(con: duckdb.DuckDBPyConnection) -> pd.DataFrame:
    """Per-species proportions (of gene_count) -- raw fractions, not standardized,
    for human-readable display (composition bars, tables). See
    get_localization_matrix for the z-scored version ordination/PERMANOVA use."""
    gene_counts = get_species_table(con)["gene_count"]
    counts = get_localization_counts(con)
    return counts.div(gene_counts, axis=0).fillna(0.0)


def get_localization_matrix(con: duckdb.DuckDBPyConnection, rows: list[str] | None = None) -> pd.DataFrame:
    """Same fractions as get_localization_fractions, standardized (z-score) per
    column for ordination/PERMANOVA (euclidean distance needs comparable scales
    across very different fraction magnitudes). `rows` subsets BEFORE z-scoring,
    same reasoning as get_assembly_matrix."""
    block = get_localization_fractions(con)
    block = _filter_rows(block, rows)
    return (block - block.mean()) / block.std(ddof=0).replace(0, 1)


def make_feature_sets(merops_families_path: Path | None = None) -> dict:
    """The registered feature sets. Built by a call (not a module constant) so the
    merops loader can carry the configured family-lookup path -- the MEROPS install
    path is per-lab/per-version, not something to hardcode at import time."""
    return {
        "assembly": {
            "loader": get_assembly_matrix, "metric": "euclidean", "method": "pcoa",
            "description": "Assembly/annotation stats (size, N50, GC%, gene count, etc.), z-scored.",
        },
        "codon": {
            "loader": get_codon_matrix, "metric": "chi-square", "method": "ca",
            "description": "Codon usage frequency profile (all 64 codons).",
        },
        "amino_acid": {
            "loader": get_aa_matrix, "metric": "chi-square", "method": "ca",
            "description": "Amino-acid usage frequency profile (all 20 residues).",
        },
        "pfam": {
            "loader": lambda con, rows=None: get_domain_matrix(con, "pfam", rows=rows, merops_families_path=merops_families_path),
            "metric": "bray-curtis", "method": "pcoa",
            "description": "Pfam domain-family content, gene-count-normalized.",
        },
        "cazy": {
            "loader": lambda con, rows=None: get_domain_matrix(con, "cazy", rows=rows, merops_families_path=merops_families_path),
            "metric": "bray-curtis", "method": "pcoa",
            "description": "CAZy carbohydrate-active-enzyme family content, gene-count-normalized.",
        },
        "merops": {
            "loader": lambda con, rows=None: get_domain_matrix(con, "merops", rows=rows, merops_families_path=merops_families_path),
            "metric": "bray-curtis", "method": "pcoa",
            "description": "MEROPS peptidase family content, gene-count-normalized.",
        },
        "localization": {
            "loader": get_localization_matrix, "metric": "euclidean", "method": "pcoa",
            "description": (
                "Composite of six protein-localization/secretion proportions, z-scored and "
                "ordinated together (not one plot per field): SignalP secretion, TMHMM "
                "transmembrane-helix count, TargetP secretion, PredGPI GPI-anchor, IDP "
                "intrinsic-disorder, WoLF PSORT extracellular."
            ),
        },
    }


def dataset_notes(n_species: int, genus_counts: dict[str, int], functional_n: int) -> str:
    """Small-n caveats text for a build, generated from the actual subset size and
    composition (never a hardcoded "11 genomes" line)."""
    from collections import Counter

    top_genera = Counter(genus_counts).most_common(5)
    genus_summary = ", ".join(f"{g} ({n})" for g, n in top_genera) if top_genera else "n/a"
    return f"""\
Dataset shape: {n_species} genomes (top genera: {genus_summary}).
`species_prefix` (in per-hit tables like pfam/cazy/merops/signalp/...) and `LOCUSTAG`
(in species/asm_stats/gene_info) are the SAME identifier for the same genome -- just
named differently by table. Always join/filter genomes through one of these two columns.

Functional-annotation coverage: {functional_n} of {n_species} genomes have hits in the
functional tables -- a genome absent from those tables will show as zeros (assembly/
BUSCO views may still be populated). Check per-table coverage before over-interpreting
an "absent family" on a genome that simply has not been screened.

Normalization choices baked into the feature matrices (do not re-derive differently):
- Pfam/CAZy/MEROPS domain-family counts are gene-count-normalized to "per 1000 genes"
  before any distance/ordination step, because raw gene counts range widely across
  genomes -- comparing raw counts would mostly measure genome/annotation size.
- A family absent from a genome is a true zero (not missing data). Where a method
  needs strictly positive values (log/fold-change), a single PSEUDOCOUNT=0.5 is added
  at matrix-build time -- never re-add a different pseudocount downstream.
- Assembly stats and the localization composite are z-scored per column (needed for
  Euclidean distance across fields with very different scales/units).
- Codon and amino-acid usage are compositional frequency profiles analyzed via
  correspondence analysis (chi-square distance), not Euclidean/Bray-Curtis.

Small-n caveats (any taxon-filtered subset is smaller still):
- PCoA/CA and PERMANOVA are unreliable at very small n: a 3-point configuration is
  exactly planar and will spuriously "explain" ~100% of variance; PERMANOVA on a
  singleton group has zero within-group variance by construction, which can look like
  a location effect when it is really a dispersion artifact. The ordination/PERMANOVA
  tools enforce minimum-n floors and attach caveat text directly in their output.
- NMDS is intentionally not the default (needs multiple random restarts to avoid a
  degenerate solution) -- offered only as an explicit opt-in, always reporting stress.
- A domain-count ordination axis can reflect assembly/annotation quality (BUSCO
  completeness, N50/contig count) or shared ancestry rather than genuine biological
  signal -- treat BUSCO/contiguity metadata and phylogenetic distance as confound
  checks, not decoration, especially for small subsets.
- `differential_families`'s "rest of clade" comparison group shrinks with any taxon
  filter -- its output reports which genomes made up "rest" for exactly this reason.
"""


# --------------------------------------------------------------------------- phylogeny
# Phyling tip labels are "{SPECIES}_{STRAIN}{suffix}" with spaces replaced by
# underscores (e.g. "Rhodotorula_mucilaginosa_DBVPG_3045_hap1.proteins"). The exact
# shape is per-dataset, so _tip_to_locustag tries several candidate spellings and
# falls back to the species name alone when STRAIN is empty.
def _tip_to_locustag(con: duckdb.DuckDBPyConnection, tip_suffix: str = ".proteins") -> dict[str, str]:
    sp = con.execute("SELECT LOCUSTAG, SPECIES, STRAIN FROM species").fetchdf()
    mapping = {}
    for row in sp.itertuples():
        base = str(row.SPECIES).replace(" ", "_")
        strain = str(row.STRAIN or "").replace(" ", "_")
        candidates = {
            f"{base}_{strain}{tip_suffix}" if strain else None,
            f"{base}{tip_suffix}",
            base,
        }
        for c in candidates:
            if c:
                mapping.setdefault(c, row.LOCUSTAG)
    return mapping


def load_tree(
    con: duckdb.DuckDBPyConnection,
    tree_path: Path | None = None,
    outgroup_genus: str | None = None,
    results_dir: Path | None = None,
    rows: list[str] | None = None,
    tip_suffix: str = ".proteins",
):
    """Bio.Phylo tree with tips renamed to LOCUSTAG, pruned to `rows` (the subset),
    and rooted on the outgroup_genus clade when one is configured.

    tree_path=None globs under results_dir via find_tree_path() rather than assuming
    a fixed filename. Tips with no matching species row are pruned rather than
    raising -- phyling builds can include extra haplotig/unphased tips (e.g. both
    hap1 and hap2 assemblies) where only one haplotype is canonical. FastTree writes
    an unrooted trifurcation at the root, so without an explicit reroot the outgroup
    placement is arbitrary -- set outgroup_genus to fix it. When `rows` is given,
    any tip outside the subset is pruned too, so a single full-cohort treefile can
    back many subset dashboards."""
    from Bio import Phylo

    if tree_path is None:
        tree_path = find_tree_path(results_dir)

    tree = Phylo.read(str(tree_path), "newick")
    tip_map = _tip_to_locustag(con, tip_suffix=tip_suffix)
    keep = set(rows) if rows else None
    unmatched = []
    for tip in tree.get_terminals():
        mapped = tip_map.get(tip.name)
        if mapped is not None and (keep is None or mapped in keep):
            tip.name = mapped
        else:
            unmatched.append(tip)
    for tip in unmatched:
        tree.prune(tip)

    if outgroup_genus:
        sp = con.execute(
            "SELECT LOCUSTAG FROM species WHERE GENUS = ?", [outgroup_genus]
        ).fetchdf()
        if sp.empty:
            raise ValueError(f"No species rows found for outgroup genus {outgroup_genus!r}")
        outgroup_locustags = set(sp["LOCUSTAG"])
        outgroup_tips = [t for t in tree.get_terminals() if t.name in outgroup_locustags]
        if not outgroup_tips:
            raise ValueError(f"Outgroup genus {outgroup_genus!r} tips not found in tree")
        tree.root_with_outgroup(outgroup_tips)
    return tree


def find_tree_path(results_dir: Path | None = None) -> Path:
    """Locate the phyling FastTree treefile under results_dir without hardcoding the
    dataset-name prefix or taxa count baked into its filename by phyling. Raises
    FileNotFoundError (caught by make_dashboard.build_payload same as any other
    tree-unavailable case) if nothing matches. If more than one lineage/build is
    present, picks the most-recently-modified match and warns -- ambiguity here
    should be visible, not silently resolved to whichever glob() returns first."""
    if results_dir is None or not results_dir.exists():
        raise FileNotFoundError(f"results_dir not found: {results_dir}")
    matches = sorted(
        results_dir.glob(TREE_GLOB), key=lambda p: p.stat().st_mtime, reverse=True
    )
    if not matches:
        raise FileNotFoundError(f"No phyling treefile found under {results_dir}/{TREE_GLOB}")
    if len(matches) > 1:
        print(
            f"  phylogeny: {len(matches)} treefiles found under {results_dir}, "
            f"using most recently modified: {matches[0]}",
            file=sys.stderr,
        )
    return matches[0]


def phylo_layout(tree) -> dict:
    """Rectangular-cladogram layout: (x, y) per node plus the segments to draw,
    precomputed in Python so the front end only needs to draw lines/circles -- no
    Newick parser in JS."""
    tree = tree.__class__(tree.root, rooted=True)  # shallow-ish; ladderize() mutates in place
    tree.ladderize()

    y_counter = {"n": 0}
    y_pos = {}
    for tip in tree.get_terminals():
        y_pos[id(tip)] = y_counter["n"]
        y_counter["n"] += 1

    def assign_y(clade):
        if clade.is_terminal():
            return y_pos[id(clade)]
        ys = [assign_y(c) for c in clade.clades]
        y = sum(ys) / len(ys)
        y_pos[id(clade)] = y
        return y

    assign_y(tree.root)

    x_pos = {id(tree.root): 0.0}

    def assign_x(clade):
        for child in clade.clades:
            x_pos[id(child)] = x_pos[id(clade)] + (child.branch_length or 0.0)
            assign_x(child)

    assign_x(tree.root)

    segments = []  # each: {x0,y0,x1,y1} -- one horizontal branch + connecting vertical per node

    def walk(clade):
        for child in clade.clades:
            segments.append(
                {
                    "x0": x_pos[id(clade)], "y0": y_pos[id(child)],
                    "x1": x_pos[id(child)], "y1": y_pos[id(child)],
                }
            )
            walk(child)
        if len(clade.clades) > 1:
            child_ys = [y_pos[id(c)] for c in clade.clades]
            segments.append(
                {
                    "x0": x_pos[id(clade)], "y0": min(child_ys),
                    "x1": x_pos[id(clade)], "y1": max(child_ys),
                }
            )

    walk(tree.root)

    tips = sorted(tree.get_terminals(), key=lambda t: y_pos[id(t)])
    return {
        "segments": segments,
        "tips": [{"locustag": t.name, "x": x_pos[id(t)], "y": y_pos[id(t)]} for t in tips],
        "tip_order": [t.name for t in tips],
        "max_x": max(x_pos.values()) if x_pos else 0.0,
    }


def gene_link(con: duckdb.DuckDBPyConnection, protein_ids: list[str]) -> pd.DataFrame:
    """Map protein_id -> gene_id/chrom/coords for 'examine this gene' drill-down links."""
    if not protein_ids:
        return pd.DataFrame(columns=["protein_id", "gene_id", "chrom", "start", "end", "strand"])
    con.execute("CREATE OR REPLACE TEMP TABLE _pids AS SELECT UNNEST(?) AS protein_id", [protein_ids])
    return con.execute(
        """
        SELECT gp.protein_id, gi.gene_id, gi.chrom, gi.start, gi.end, gi.strand
        FROM _pids p
        JOIN gene_proteins gp ON gp.protein_id = p.protein_id
        JOIN gene_info gi ON gi.gene_id = gp.gene_id
        """
    ).fetchdf()
