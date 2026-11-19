"""DAP-seq target sets of Lax et al. 2025 for R. microsporus ATCC 11559, mapped to BFD (FBF4BF29), and their
conservation across the genus-representative proteomes of BFD.

Pipeline (see workflows/tf_targets.sh): prep -> map -> repdb -> conserve -> tables -> figures. Each step skips
when its output exists. The pure functions (target_sets, count_taxa, compare_to_background, profile_targets)
are unit tested with toy tables.

Mapping rate is not a conservation measure. A JGI target without a BFD partner may come from a different gene
model of another assembly, so it is not evidence of a lineage-specific gene.
"""
import gzip
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as st

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "summary"))
import jgi_map as jm  # noqa: E402
from stats import bh_adjust  # noqa: E402

BFD_GENOME = "FBF4BF29"
MIN_TARGETS = 20
HI_ORDERS = 100          # "widely conserved" cut for the fraction columns and the figure
N_BINS = 5               # protein-length strata
MIN_STRATUM = 5          # targets needed in a stratum to report its own effect
MAX_SEQS = 3000
MIN_COV = 0.5
EVALUE = 1e-5
RUNS = Path("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs")
RES = Path("/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/resources")
JGI = RES / "jgi_atcc_11559/Rhimi59_2/Mycocosm/Annotation"
JGI_FASTA = JGI / "All_models,_Filtered_and_Not/Proteins/Rhimi59_2_all_proteins_20220409.aa.fasta.gz"
JGI_GFF = JGI / "Filtered_Models___best__/Genes/Rhimi59_2_GeneCatalog_20220409.gff3.gz"
MOESM4 = RES / "s41467-025-65177-1_suppl/41467_2025_65177_MOESM4_ESM.xlsx"
MOESM6 = RES / "s41467-025-65177-1_suppl/41467_2025_65177_MOESM6_ESM.xlsx"


# ---------------------------------------------------------------- target sets
def target_sets(targets: pd.DataFrame, id_map: pd.DataFrame) -> pd.DataFrame:
    """Targets with their BFD protein. Unmapped targets are kept with an empty bfd_protein_id."""
    m = id_map[["jgi_protein_id", "bfd_protein_id"]].copy()
    m["jgi_protein_id"] = m["jgi_protein_id"].astype(str)
    m = m[m["bfd_protein_id"].fillna("") != ""].drop_duplicates("jgi_protein_id")
    t = targets[["tf_name", "jgi_protein_id"]].copy()
    t["jgi_protein_id"] = t["jgi_protein_id"].astype(str)
    out = t.merge(m, on="jgi_protein_id", how="left")
    out["bfd_protein_id"] = out["bfd_protein_id"].fillna("")
    return out[["tf_name", "jgi_protein_id", "bfd_protein_id"]]


def set_sizes(ts: pd.DataFrame) -> pd.DataFrame:
    """Per TF: number of targets, mapped and unmapped."""
    g = ts.assign(_m=ts["bfd_protein_id"].fillna("") != "").groupby("tf_name", sort=False)["_m"]
    out = pd.DataFrame({"n_targets": g.size(), "n_mapped": g.sum().astype(int)}).reset_index()
    out["n_unmapped"] = out["n_targets"] - out["n_mapped"]
    return out


# ---------------------------------------------------------------- conservation
def count_taxa(hits: pd.DataFrame, taxonomy: pd.DataFrame, queries, min_cov=MIN_COV, evalue=EVALUE) -> pd.DataFrame:
    """Distinct phyla, orders and genera among the hits that pass the e-value and query coverage thresholds.
    The genome of a hit is split_part(target, '_', 1). Genomes whose order (or phylum, genus) is empty add
    nothing to that count. Queries without a kept hit get 0."""
    h = hits[(hits["evalue"] <= evalue) & (hits["qcov"] >= min_cov)].copy()
    h["species_prefix"] = h["target"].astype(str).str.split("_", n=1).str[0]
    tax = taxonomy[["species_prefix", "PHYLUM", "ORDER", "GENUS"]].drop_duplicates("species_prefix")
    h = h.merge(tax, on="species_prefix", how="inner")
    g = h.groupby("query")
    cnt = pd.DataFrame({"n_phyla": g["PHYLUM"].nunique(), "n_orders": g["ORDER"].nunique(), "n_genera": g["GENUS"].nunique()})
    out = pd.DataFrame({"bfd_protein_id": list(queries)}).merge(cnt, left_on="bfd_protein_id", right_index=True, how="left")
    for c in ("n_phyla", "n_orders", "n_genera"):
        out[c] = out[c].fillna(0).astype(int)
    out.attrs.update({"min_cov": min_cov, "evalue": evalue})
    return out


def extract_fasta(src, ids, dst):
    """Write the records of src (plain or .gz) whose first header word is in ids. Header is cut to that word."""
    ids = set(ids)
    n, keep = 0, False
    with jm._open_text(src) as f, open(dst, "w") as o:
        for line in f:
            if line.startswith(">"):
                pid = line[1:].split()[0]
                keep = pid in ids
                if keep:
                    o.write(f">{pid}\n")
                    n += 1
            elif keep:
                o.write(line if line.endswith("\n") else line + "\n")
    return n


def search_genus_reps(query_fasta, db, hits_tsv, tmp, min_cov=MIN_COV, evalue=EVALUE, threads=None,
                      split_memory_limit="30G", max_seqs=MAX_SEQS, sensitivity=4.0):
    """MMseqs2 search of query_fasta against the genus-representative database. Query coverage is
    --cov-mode 2. max_seqs caps the hits per query; sensitivity is below the default 5.7 to fit the time limit."""
    exe = jm.need_mmseqs()
    tmp = Path(tmp)
    tmp.mkdir(parents=True, exist_ok=True)
    hits_tsv = Path(hits_tsv)
    part = hits_tsv.with_name(hits_tsv.name + ".part")     # a killed run leaves .part, never a complete-looking file
    cmd = [exe, "easy-search", str(query_fasta), str(db), str(part), str(tmp / "mm_tmp"),
           "-e", str(evalue), "-c", str(min_cov), "--cov-mode", "2", "-s", str(sensitivity),
           "--max-seqs", str(max_seqs), "--format-output", "query,target,evalue,qcov,bits",
           "--threads", str(threads or jm.threads_default()), "--split-memory-limit", str(split_memory_limit), "-v", "1"]
    subprocess.run(cmd, check=True)
    part.replace(hits_tsv)


def read_conservation_hits(path):
    return pd.read_csv(path, sep="\t", header=None, names=["query", "target", "evalue", "qcov", "bits"],
                       dtype={"query": str, "target": str})


def conservation(target_bfd_ids, genus_reps_proteome: Path, out: Path, scratch: Path, min_cov=MIN_COV, evalue=EVALUE,
                 *, query_fasta=None, taxonomy=None, hits=None, threads=None, hits_tsv=None) -> pd.DataFrame:
    """Per BFD protein: n_phyla, n_orders, n_genera among the genus-representative proteomes with a hit
    (e-value <= evalue, query coverage >= min_cov). genus_reps_proteome is the MMseqs2 database (or FASTA).
    hits: a hit table (DataFrame or path) to use instead of a search. Thresholds are on the first line of out."""
    ids = list(dict.fromkeys(target_bfd_ids))
    if hits is None:
        if query_fasta is None:
            raise ValueError("query_fasta is needed when no hit table is given")
        scratch = Path(scratch)
        scratch.mkdir(parents=True, exist_ok=True)
        qf = scratch / "conservation_queries.faa"
        n = extract_fasta(query_fasta, ids, qf)
        if n != len(ids):
            raise ValueError(f"{len(ids) - n} query proteins are not in {query_fasta}")
        hits = Path(hits_tsv) if hits_tsv else scratch / "conservation_hits.tsv"
        if not (hits.exists() and hits.stat().st_size > 0):
            search_genus_reps(qf, genus_reps_proteome, hits, scratch, min_cov, evalue, threads)
    if not isinstance(hits, pd.DataFrame):
        hits = read_conservation_hits(hits)
    df = count_taxa(hits, taxonomy, ids, min_cov, evalue)
    out = Path(out)
    with open(out, "w") as fh:
        fh.write(f"# min_cov={min_cov} evalue={evalue} cov_mode=2 (query) taxa=distinct PHYLUM/ORDER/GENUS of genus representatives\n")
        df.to_csv(fh, sep="\t", index=False)
    return df


# ---------------------------------------------------------------- comparison to background
def compare_to_background(conservation_df: pd.DataFrame, tf_targets: pd.DataFrame, background_ids) -> pd.DataFrame:
    """Per TF: n_orders of its mapped targets against n_orders of all background proteins (two-sided
    Mann-Whitney U). effect is the rank-biserial correlation 2U/(n1*n2)-1, positive when targets are more widely
    conserved. BH over the tested TFs. TFs with fewer than 20 mapped targets are not tested (reason says so)."""
    cons = conservation_df.set_index("bfd_protein_id")["n_orders"]
    bg = cons.reindex(list(background_ids)).dropna().to_numpy(float)
    rows = []
    for tf, g in tf_targets.groupby("tf_name", sort=False):
        mapped = sorted(set(g.loc[g["bfd_protein_id"].fillna("") != "", "bfd_protein_id"]))
        vals = cons.reindex(mapped).dropna().to_numpy(float)
        row = {"tf_name": tf, "n_targets_bfd": len(mapped), "median_n_orders": float(np.median(vals)) if len(vals) else np.nan,
               "background_median_n_orders": float(np.median(bg)), "effect": np.nan, "p": np.nan, "reason": "",
               "frac_n_orders_ge_100": float(np.mean(vals >= HI_ORDERS)) if len(vals) else np.nan,
               "background_frac_n_orders_ge_100": float(np.mean(bg >= HI_ORDERS))}
        if len(vals) < MIN_TARGETS:
            row["reason"] = f"fewer than {MIN_TARGETS} mapped targets"
        else:
            u, p = st.mannwhitneyu(vals, bg, alternative="two-sided")
            row["effect"] = float(2 * u / (len(vals) * len(bg)) - 1)
            row["p"] = float(p)
        rows.append(row)
    out = pd.DataFrame(rows)
    out["p_adj"] = bh_adjust(out["p"].to_numpy())
    return out[["tf_name", "n_targets_bfd", "median_n_orders", "background_median_n_orders", "effect", "p", "p_adj", "reason",
                "frac_n_orders_ge_100", "background_frac_n_orders_ge_100"]]


def length_bins(lengths, bg_ids, n_bins=N_BINS):
    """Stratum 1..n_bins of each protein. Edges are the quantiles of the background (all FBF4BF29) lengths."""
    bg_len = np.asarray([lengths[i] for i in bg_ids], float)
    edges = np.quantile(bg_len, np.arange(1, n_bins) / n_bins)
    return lambda ids: np.searchsorted(edges, np.asarray([lengths[i] for i in ids], float), side="right") + 1, edges


def _van_elteren(groups):
    """Stratified Wilcoxon rank-sum (van Elteren). groups: list of (targets, background) value arrays per stratum.
    Stratum weight 1/(N+1). Returns (z, two-sided p); z > 0 when targets rank higher. Tie-corrected variance."""
    num = var = 0.0
    for x, y in groups:
        n1, n2 = len(x), len(y)
        if n1 == 0 or n2 == 0:
            continue
        N = n1 + n2
        allv = np.concatenate([x, y])
        ranks = st.rankdata(allv)
        w = 1.0 / (N + 1)
        _, cnt = np.unique(allv, return_counts=True)
        tie = float(np.sum(cnt ** 3 - cnt)) / (N * (N - 1)) if N > 1 else 0.0
        v = n1 * n2 / 12.0 * ((N + 1) - tie)
        num += w * (ranks[:n1].sum() - n1 * (N + 1) / 2.0)
        var += w * w * v
    if var <= 0:
        return np.nan, np.nan
    z = num / np.sqrt(var)
    return float(z), float(2 * st.norm.sf(abs(z)))


def _weighted_median(values, weights):
    o = np.argsort(values, kind="stable")
    v, w = np.asarray(values, float)[o], np.asarray(weights, float)[o]
    c = np.cumsum(w)
    return float(v[np.searchsorted(c, 0.5 * c[-1])])


def stratified_tests(conservation_df: pd.DataFrame, tf_targets: pd.DataFrame, background_ids, lengths, n_bins=N_BINS):
    """Length-stratified comparison of each TF's mapped targets with the background (all FBF4BF29 proteins).
    Strata: quintiles of the background protein length. Per stratum: rank-biserial effect (2U/(n1 n2) - 1), reported
    when the stratum holds >= 5 targets. Pooled effect_strat = mean of the per-stratum effects over strata with >= 5
    targets, weighted by the number of targets in the stratum. p_strat = van Elteren stratified rank-sum test (weights
    1/(N+1)), two-sided. BH over the tested TFs (>= 20 mapped targets). expected_median_n_orders: weighted median of
    the background where each background protein is weighted by (targets in its stratum / background in its stratum),
    that is the median expected if the targets had the background's n_orders at their own length profile.
    Returns (tests, strata): tests has one row per TF; strata has one row per TF and stratum."""
    cons = conservation_df.set_index("bfd_protein_id")["n_orders"]
    bg_ids = [i for i in dict.fromkeys(background_ids) if i in cons.index and i in lengths.index]
    bin_of, edges = length_bins(lengths, bg_ids, n_bins)
    bg_vals = cons.reindex(bg_ids).to_numpy(float)
    bg_bin = bin_of(bg_ids)
    rows, srows = [], []
    for tf, g in tf_targets.groupby("tf_name", sort=False):
        mapped = [i for i in sorted(set(g.loc[g["bfd_protein_id"].fillna("") != "", "bfd_protein_id"])) if i in cons.index and i in lengths.index]
        vals = cons.reindex(mapped).to_numpy(float)
        tb = bin_of(mapped) if mapped else np.array([], int)
        row = {"tf_name": tf, "n_targets_bfd": len(mapped), "effect_strat": np.nan, "z_strat": np.nan, "p_strat": np.nan,
               "expected_median_n_orders": np.nan, "reason": ""}
        groups, effs, wts = [], [], []
        for k in range(1, n_bins + 1):
            x, y = vals[tb == k], bg_vals[bg_bin == k]
            eff = np.nan
            if len(x) >= MIN_STRATUM and len(y):
                u, _ = st.mannwhitneyu(x, y, alternative="two-sided")
                eff = float(2 * u / (len(x) * len(y)) - 1)
                effs.append(eff)
                wts.append(len(x))
            groups.append((x, y))
            srows.append({"tf_name": tf, "stratum": k, "length_min": float(edges[k - 2]) if k > 1 else np.nan,
                          "length_max": float(edges[k - 1]) if k < n_bins else np.nan, "n_targets": int(len(x)),
                          "n_background": int(len(y)), "median_targets": float(np.median(x)) if len(x) else np.nan,
                          "median_background": float(np.median(y)) if len(y) else np.nan, "effect": eff})
        if len(mapped) < MIN_TARGETS:
            row["reason"] = f"fewer than {MIN_TARGETS} mapped targets"
        else:
            if effs:
                row["effect_strat"] = float(np.average(effs, weights=wts))
            row["z_strat"], row["p_strat"] = _van_elteren(groups)
            nt = np.bincount(tb, minlength=n_bins + 1)[1:]
            nb = np.bincount(bg_bin, minlength=n_bins + 1)[1:]
            w = (nt / np.maximum(nb, 1))[bg_bin - 1]
            row["expected_median_n_orders"] = _weighted_median(bg_vals, w)
        rows.append(row)
    out = pd.DataFrame(rows)
    out["p_adj_strat"] = bh_adjust(out["p_strat"].to_numpy())
    cols = ["tf_name", "n_targets_bfd", "effect_strat", "z_strat", "p_strat", "p_adj_strat", "expected_median_n_orders", "reason"]
    return out[cols], pd.DataFrame(srows)


def fasta_lengths(path) -> pd.Series:
    """Protein ID -> length from a FASTA (plain or .gz)."""
    d, pid, n = {}, None, 0
    with jm._open_text(path) as f:
        for line in f:
            if line.startswith(">"):
                if pid is not None:
                    d[pid] = n
                pid, n = line[1:].split()[0], 0
            else:
                n += len(line.strip())
    if pid is not None:
        d[pid] = n
    return pd.Series(d, dtype=float)


# ---------------------------------------------------------------- profile
def profile_targets(target_sets: pd.DataFrame, pfam_hits: pd.DataFrame, swissprot: pd.DataFrame, top=5) -> pd.DataFrame:
    """Per TF: top Pfam families and top SwissProt best hits among the mapped targets, by number of proteins.
    pfam_hits: protein_id, pfam_id, pfam_acc (version stripped here). swissprot: protein_id, accession,
    entry_name, protein_name. Ties are broken by term."""
    ts = target_sets[target_sets["bfd_protein_id"].fillna("") != ""][["tf_name", "bfd_protein_id"]].drop_duplicates()
    n_t = ts.groupby("tf_name").size()
    pf = pfam_hits.copy()
    pf["term"] = pf["pfam_acc"].astype(str).str.split(".").str[0]
    pf["term_name"] = pf["pfam_id"]
    sp = swissprot.copy()
    sp["term"] = sp["accession"]
    sp["term_name"] = sp["entry_name"].astype(str) + " " + sp["protein_name"].astype(str)
    parts = []
    for source, tab in (("pfam", pf), ("swissprot", sp)):
        j = ts.merge(tab[["protein_id", "term", "term_name"]], left_on="bfd_protein_id", right_on="protein_id")
        c = j.groupby(["tf_name", "term"], sort=False).agg(term_name=("term_name", "first"),
                                                           n_proteins=("bfd_protein_id", "nunique")).reset_index()
        c = c.sort_values(["tf_name", "n_proteins", "term"], ascending=[True, False, True])
        c["rank"] = c.groupby("tf_name").cumcount() + 1
        c = c[c["rank"] <= top].copy()
        c.insert(1, "source", source)
        parts.append(c)
    out = pd.concat(parts, ignore_index=True)
    out["n_targets_bfd"] = out["tf_name"].map(n_t).astype(int)
    out["frac"] = out["n_proteins"] / out["n_targets_bfd"]
    return out[["tf_name", "source", "rank", "term", "term_name", "n_proteins", "n_targets_bfd", "frac"]].sort_values(
        ["tf_name", "source", "rank"]).reset_index(drop=True)


def load_profile_sources(tables_dir, ids, scratch, min_neglog_e=5.0):
    """Pfam rows and SwissProt best hits of the given BFD protein IDs, read with DuckDB through a join against a
    small ID table. SwissProt hits need neg_log10_evalue >= min_neglog_e (e-value <= 1e-5)."""
    import duckdb
    tables_dir = Path(tables_dir)
    con = duckdb.connect()
    con.execute(f"SET temp_directory='{scratch}'")
    con.register("ids_df", pd.DataFrame({"protein_id": sorted(set(ids))}))
    con.execute("CREATE TEMP TABLE ids AS SELECT * FROM ids_df")
    pf = con.execute(f"SELECT DISTINCT p.protein_id, p.pfam_id, p.pfam_acc FROM read_parquet('{tables_dir}/pfam.parquet') p "
                     "JOIN ids USING (protein_id)").df()
    sp = con.execute(f"SELECT s.protein_id, s.accession, s.entry_name, s.protein_name FROM "
                     f"read_parquet('{tables_dir}/swissprot_best.parquet') s JOIN ids USING (protein_id) "
                     f"WHERE s.neg_log10_evalue >= {min_neglog_e}").df()
    return pf, sp


# ---------------------------------------------------------------- genus representatives
def genus_reps(summary_dir, suppress=None):
    """Genus representatives by the ordination rule (pick_complete_reps: best BUSCO, then N50, then lowest id,
    from genomes with data in the codon, amino-acid and gene-structure tables), after the suppress list.
    Returns (reps, genome_all after suppression)."""
    from suppress import drop_suppressed, load_suppress
    from summary_ordination import pick_complete_reps
    sd = Path(summary_dir)
    sup = suppress if isinstance(suppress, pd.DataFrame) else load_suppress(suppress)

    def rd(n):
        return drop_suppressed(pd.read_parquet(sd / f"{n}.parquet"), sup)[0]
    ga = rd("genome_all")
    mats = {"codon": rd("codon_genome"), "aa": rd("aa_genome"), "structure": rd("gene_structure_genome")}
    return pick_complete_reps(ga, mats), ga


def proteome_paths(genomes_tsv):
    """species_prefix -> protein FASTA path from seqdb/genomes.tsv (third column)."""
    d = {}
    with open(genomes_tsv) as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) >= 3:
                d[c[0]] = c[2]
    return d


# ---------------------------------------------------------------- figures
CAVEAT = ("Hypothesis only. n_orders counts orders with a homology hit (domains, paralogs), not orthologs. Targets overlap "
          "across TFs, so the tests are not independent. The background (all genes) contains the targets.")


def family_colors(families):
    import matplotlib.pyplot as plt
    pal = list(plt.get_cmap("tab20").colors) + list(plt.get_cmap("tab20b").colors)
    return {f: pal[i % len(pal)] for i, f in enumerate(sorted(set(families)))}


def _fig_conservation(tests, ecdf):
    """Top: ECDF of n_orders (pooled mapped targets, pooled without Znf1, background). Bottom left: fraction of each TF's
    targets with n_orders >= 100, dashed line = background fraction. Bottom right: length-stratified effect (filled when
    BH p < 0.05 on the stratified test). n_orders is bimodal (many proteins with ~0 or ~190 orders), so no boxes."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    ok = tests[tests["reason"].fillna("") == ""].sort_values("frac_n_orders_ge_100").reset_index(drop=True)
    if ok.empty:
        return None
    cm = family_colors(ok["tf_family"].fillna("unknown"))
    h = max(7.0, 0.19 * len(ok) + 4.6)
    fig = plt.figure(figsize=(9.5, h))
    gs = fig.add_gridspec(2, 2, height_ratios=[2.2, 0.19 * len(ok) + 1.0], hspace=0.28, wspace=0.06, left=0.17, right=0.78, top=0.96, bottom=0.07)
    ae = fig.add_subplot(gs[0, :])
    for col, lab, c, ls in (("ecdf_background", "all FBF4BF29 proteins", "#c0392b", "--"), ("ecdf_targets", "mapped targets (pooled, unique proteins)", "#1f4e79", "-"),
                            ("ecdf_targets_no_Znf1", "mapped targets without Znf1", "#6fa8dc", "-")):
        ae.step(ecdf["n_orders"], ecdf[col], where="post", color=c, ls=ls, lw=1.6, label=lab)
    ae.axvline(HI_ORDERS, color="grey", lw=0.6, ls=":")
    ae.set_xlabel("n_orders (orders with a hit among genus representatives)", fontsize=8)
    ae.set_ylabel("cumulative fraction", fontsize=8)
    ae.legend(fontsize=7, frameon=False, loc="upper left")
    ae.set_title("Conservation of DAP-seq targets (TFs with >= 20 mapped targets)", fontsize=9)
    ab = fig.add_subplot(gs[1, 0])
    ac = fig.add_subplot(gs[1, 1], sharey=ab)
    y = np.arange(len(ok))
    cols = [cm[f] for f in ok["tf_family"].fillna("unknown")]
    ab.scatter(ok["frac_n_orders_ge_100"], y, c=cols, s=18, zorder=3)
    ab.axvline(float(ok["background_frac_n_orders_ge_100"].iloc[0]), color="#c0392b", ls="--", lw=1.2)
    ab.set_yticks(y)
    ab.set_yticklabels([f"{t} (n={n})" for t, n in zip(ok["tf_name"], ok["n_targets_bfd"])], fontsize=6)
    ab.set_xlabel(f"fraction of targets with n_orders >= {HI_ORDERS}", fontsize=8)
    ab.grid(axis="x", lw=0.3)
    sig = (ok["p_adj_strat"] < 0.05).to_numpy()
    ac.scatter(ok["effect_strat"][sig], y[sig], c=[c for c, s_ in zip(cols, sig) if s_], s=18, zorder=3)
    ac.scatter(ok["effect_strat"][~sig], y[~sig], facecolors="none", edgecolors=[c for c, s_ in zip(cols, sig) if not s_], s=18, zorder=3)
    ac.axvline(0, color="grey", lw=0.8)
    ac.set_xlabel("length-stratified effect (rank-biserial)", fontsize=8)
    ac.grid(axis="x", lw=0.3)
    plt.setp(ac.get_yticklabels(), visible=False)
    handles = [Line2D([0], [0], marker="o", ls="", color=c, label=f, markersize=5) for f, c in cm.items()]
    handles += [Line2D([0], [0], marker="o", ls="", mfc="none", mec="k", label="open: stratified BH p >= 0.05", markersize=5),
                Line2D([0], [0], marker="o", ls="", color="k", label="filled: stratified BH p < 0.05", markersize=5)]
    fig.legend(handles=handles, fontsize=6, title="paper TF family", title_fontsize=7, loc="center left", bbox_to_anchor=(0.79, 0.4), frameon=False)
    fig.text(0.01, 0.005, CAVEAT, fontsize=5.8, wrap=True, va="bottom")
    return fig


def _fig_sizes(sizes):
    import matplotlib.pyplot as plt
    s = sizes.sort_values("n_targets")
    fig, ax = plt.subplots(figsize=(7.5, max(4.0, 0.15 * len(s) + 1.6)))
    y = np.arange(len(s))
    ax.barh(y, s["n_mapped"], color="#2a6f97", label="mapped to BFD (RBH, best models)")
    ax.barh(y, s["n_unmapped"], left=s["n_mapped"], color="#c9c9c9", label="no BFD partner")
    ax.set_yticks(y)
    ax.set_ylim(-0.7, len(s) - 0.3)
    ax.set_yticklabels(s["tf_name"], fontsize=5)
    ax.set_xlabel("DAP-seq targets (JGI protein IDs)")
    ax.set_title("Targets with and without a BFD partner, per TF", fontsize=9)
    ax.legend(fontsize=7, frameon=False, loc="lower right")
    fig.tight_layout()
    return fig


def read_commented(path, **kw):
    return pd.read_csv(path, sep="\t", comment="#", keep_default_na=False, na_values=[""], **kw)


def page_tf_targets(out_dir, tables_dir=None, o=None):
    """Draw T21_target_conservation and T21_target_set_sizes from the T21 tables (T21_target_background_tests,
    T21_target_ecdf, T21_tf_set_sizes). tables_dir default: <out_dir>/tables_tf. Figures go to <out_dir>/figures/tf
    (Out level 'tf'). Returns the list of Made records."""
    from figs_common import safe
    from make_figures import Out
    out_dir = Path(out_dir)
    td = Path(tables_dir) if tables_dir else out_dir / "tables_tf"
    o = o or Out(out_dir, "tf")
    tests = read_commented(td / "T21_target_background_tests.tsv")
    tests["reason"] = tests["reason"].fillna("")
    ecdf = read_commented(td / "T21_target_ecdf.tsv")
    sizes = read_commented(td / "T21_tf_set_sizes.tsv")
    ntest = int((tests["reason"] == "").sum())
    o.info = {"n_genomes": None, "n_groups": ntest, "pooled_share": 0.0}
    o.meta["T21_target_conservation"] = {"signal": ntest, "signal_kind": "n_tfs_tested"}
    safe(o, "T21_target_conservation", "tf", lambda: _fig_conservation(tests, ecdf))
    safe(o, "T21_target_set_sizes", "tf", lambda: _fig_sizes(sizes))
    return o.made


# ---------------------------------------------------------------- pipeline steps
def _gz_write(df, path):
    part = Path(str(path) + ".part")
    with gzip.open(part, "wt") as f:
        df.to_csv(f, sep="\t", index=False)
    part.replace(path)


def _write_fasta_tmp(path: Path, records):
    tmp = path.with_name(path.name + ".tmp")
    opener = gzip.open if path.suffix == ".gz" else open
    n = 0
    with opener(tmp, "wt") as o:
        for pid, seq in records:
            o.write(f">{pid}\n{seq}\n")
            n += 1
    tmp.replace(path)
    return n


def step_prep(stage: Path, bfd_src=None, jgi_fasta=None, jgi_gff=None):
    """Write the FBF4BF29 proteome, the best-model JGI proteins (the proteinIds of the filtered-models GFF3) and the
    JGI proteins of all models (gz), all with plain ID headers. Each file is skipped when it exists."""
    stage.mkdir(parents=True, exist_ok=True)
    bfd = stage / "FBF4BF29.faa"
    if not bfd.exists():
        src = bfd_src or proteome_paths(RUNS / "seqdb/genomes.tsv")[BFD_GENOME]
        n = 0
        tmp = bfd.with_suffix(".tmp")
        with jm._open_text(src) as f, open(tmp, "w") as o:
            for line in f:
                if line.startswith(">"):
                    o.write(f">{line[1:].split()[0]}\n")
                    n += 1
                else:
                    o.write(line)
        tmp.replace(bfd)
        print(f"wrote {bfd} ({n} proteins) from {src}", flush=True)
    allf, bestf = stage / "Rhimi59_2_all_models.faa.gz", stage / "Rhimi59_2_best_models.faa"
    if not (allf.exists() and bestf.exists()):
        prot = jm.read_jgi_proteins(jgi_fasta or JGI_FASTA)
        best = list(jm.read_jgi_catalog(jgi_gff or JGI_GFF)["protein_id"])
        miss = [b for b in best if b not in prot]
        if miss:
            raise ValueError(f"{len(miss)} best-model protein IDs are not in the protein FASTA, e.g. {miss[:5]}")
        if not allf.exists():
            print(f"wrote {allf} ({_write_fasta_tmp(allf, ((k, v[1]) for k, v in prot.items()))} proteins)", flush=True)
        if not bestf.exists():
            print(f"wrote {bestf} ({_write_fasta_tmp(bestf, ((k, prot[k][1]) for k in best))} proteins)", flush=True)


# tag -> (which JGI set, min_pident, min_cov). best = filtered best models on both sides (main). all = all 258,092 models.
MAPS = {"best_p90_c0.8": ("best", 90.0, 0.8),      # main
        "best_p70_c0.7": ("best", 70.0, 0.7),      # sensitivity
        "best_p90_c0.5": ("best", 90.0, 0.5),      # sensitivity: identity kept, coverage relaxed
        "all_p90_c0.8": ("all", 90.0, 0.8)}        # secondary: reverse search against all models
MAIN_MAP = "best_p90_c0.8"


def jgi_fasta_for(stage: Path, which):
    return stage / ("Rhimi59_2_best_models.faa" if which == "best" else "Rhimi59_2_all_models.faa.gz")


def step_map(stage: Path, scratch: Path):
    """RBH maps in MAPS. Raw forward and reverse hit tables are kept in <stage>/raw_hits as gz."""
    for tag, (which, pid, cov) in MAPS.items():
        out = stage / f"jgi_to_bfd_{tag}.tsv"
        if out.exists():
            continue
        jm.map_jgi_to_bfd(jgi_fasta_for(stage, which), stage / "FBF4BF29.faa", out, scratch / f"map_{tag}", pid, cov,
                          raw_dir=stage / "raw_hits", tag=tag)
        print(f"wrote {out}", flush=True)


def step_repdb(stage: Path, summary_dir):
    """Concatenate the genus-representative proteomes into one MMseqs2 database (built once)."""
    db = stage / "genus_reps_db" / "genus_reps"
    info = stage / "genus_reps_db.info.tsv"
    if info.exists() and (db.parent / "genus_reps.dbtype").exists():
        return db
    reps, _ = genus_reps(summary_dir)
    paths = proteome_paths(RUNS / "seqdb/genomes.tsv")
    db.parent.mkdir(parents=True, exist_ok=True)
    exe = jm.need_mmseqs()
    miss = [g for g in reps["species_prefix"] if g not in paths or not Path(paths[g]).exists()]
    if miss:
        print(f"WARNING {len(miss)} representatives without a protein file: {miss[:10]}", flush=True)
    use = [g for g in reps["species_prefix"] if g not in set(miss)]
    reps.loc[reps["species_prefix"].isin(use), ["species_prefix", "PHYLUM", "ORDER", "GENUS"]].to_csv(
        stage / "genus_reps.tsv", sep="\t", index=False)
    p = subprocess.Popen([exe, "createdb", "stdin", str(db), "--dbtype", "1", "-v", "1"], stdin=subprocess.PIPE)
    for g in use:
        with jm._open_text(paths[g]) as f:
            for line in f:
                p.stdin.write((line if line.endswith("\n") else line + "\n").encode())
    p.stdin.close()
    if p.wait() != 0:
        raise RuntimeError("mmseqs createdb failed")
    n = sum(1 for _ in open(str(db) + ".index"))
    size = sum(x.stat().st_size for x in db.parent.iterdir())
    pd.DataFrame([{"n_genomes": len(use), "n_missing": len(miss), "n_proteins": n, "db_bytes": size}]).to_csv(
        info, sep="\t", index=False)
    print(f"genus reps DB: {len(use)} genomes, {n} proteins, {size / 1e9:.1f} GB", flush=True)
    return db


def _taxonomy(stage, summary_dir):
    _, ga = genus_reps(summary_dir)
    return ga


def step_conserve(stage: Path, scratch: Path, summary_dir):
    """Search ALL FBF4BF29 proteins against the genus-representative DB; count taxa per protein."""
    out = stage / "T21_conservation_all_bfd.tsv"
    if out.exists():
        return
    bfd = stage / "FBF4BF29.faa"
    ids = jm.fasta_ids(bfd)
    hits_gz = stage / "conservation_hits.tsv.gz"
    hits = hits_gz if hits_gz.exists() else stage / "conservation_hits.tsv"
    ga = _taxonomy(stage, summary_dir)
    conservation(ids, stage / "genus_reps_db" / "genus_reps", out.with_suffix(".tmp"), scratch, query_fasta=bfd,
                 taxonomy=ga, hits_tsv=hits)
    out.with_suffix(".tmp").rename(out)
    if not hits_gz.exists() and hits.exists():
        subprocess.run(["gzip", "-f", str(hits)], check=True)


def step_tables(stage: Path, tables_tf: Path, scratch: Path, runs_tables: Path):
    tables_tf.mkdir(parents=True, exist_ok=True)
    names = jm.read_tf_names(HERE.parent.parent / "config" / "paper_tf_names.tsv")
    targets = jm.read_targets(MOESM4, names)
    tested = jm.read_tested_tfs(MOESM6)
    best = set(jm.read_jgi_catalog(JGI_GFF)["protein_id"])
    maps = {tag: jm.read_map(stage / f"jgi_to_bfd_{tag}.tsv") for tag in MAPS}
    uniq_targets, tested_ids = targets["jgi_protein_id"].unique(), tested["jgi_protein_id"].unique()
    # ---- mapping rates: main map first; the all-models rows are labelled and never given a one-to-one rate
    rates = []
    for tag, (which, pid, cov) in MAPS.items():
        m = maps[tag]
        lab = "best_models_both_sides" if which == "best" else "all_models_reverse"
        for what, ids in (("tested_tfs", tested_ids), ("targets", uniq_targets)):
            r = jm.mapping_rate(m, ids, best_ids=best)
            rates.append({"map": lab, "thresholds": f"p{pid:g}_c{cov:g}", "set": what, "n_ids": r["n_ids"], "n_mapped": r["n_mapped"],
                          "rate": r["rate"], "n_distinct_bfd": r["n_distinct_bfd"], "main": tag == MAIN_MAP,
                          "note": "" if which == "best" else "reverse search against all models: an alternate model can break reciprocity"})
        r = jm.mapping_rate(m, list(m["jgi_protein_id"]))
        rates.append({"map": lab, "thresholds": f"p{pid:g}_c{cov:g}", "set": "all_models_in_fasta" if which == "all" else "best_models_all",
                      "n_ids": r["n_ids"], "n_mapped": r["n_mapped"], "rate": float("nan") if which == "all" else r["rate"],
                      "n_distinct_bfd": r["n_distinct_bfd"], "main": False,
                      "note": ("not a one-to-one rate: many JGI models map to the same BFD protein; use n_distinct_bfd" if which == "all" else "")})
    pd.DataFrame(rates).to_csv(tables_tf / "T21_jgi_map_rate.tsv", sep="\t", index=False)
    # ---- target sets (main map) and set sizes
    ts = target_sets(targets, maps[MAIN_MAP])
    tinfo = targets.drop_duplicates("tf_name")[["tf_name", "tf_family", "tf_jgi_id"]].merge(
        names[["jgi_id", "success_frip5_any"]], left_on="tf_jgi_id", right_on="jgi_id", how="left").drop(columns=["jgi_id"])
    ts = ts.merge(tinfo[["tf_name", "tf_family"]], on="tf_name", how="left")
    ts["in_best_models"] = ts["jgi_protein_id"].isin(best)
    _gz_write(ts, tables_tf / "T21_target_sets.tsv.gz")
    ss = set_sizes(ts)
    ss["n_bfd_distinct"] = ss["tf_name"].map(ts[ts["bfd_protein_id"] != ""].groupby("tf_name")["bfd_protein_id"].nunique()).fillna(0).astype(int)
    for tag in ("best_p70_c0.7", "best_p90_c0.5"):
        t2 = target_sets(targets, maps[tag])
        ss[f"n_mapped_{tag}"] = ss["tf_name"].map(t2[t2["bfd_protein_id"] != ""].groupby("tf_name").size()).fillna(0).astype(int)
    ss = ss.merge(tinfo, on="tf_name", how="left")
    ss.to_csv(tables_tf / "T21_tf_set_sizes.tsv", sep="\t", index=False)
    # ---- conservation of the targets and the background
    allc = pd.read_csv(stage / "T21_conservation_all_bfd.tsv", sep="\t", comment="#")
    bg = jm.fasta_ids(stage / "FBF4BF29.faa")
    lengths = fasta_lengths(stage / "FBF4BF29.faa")
    mapped = ts[ts["bfd_protein_id"] != ""]
    capped = _capped_queries(stage)
    st_tests, strat = stratified_tests(allc, ts, bg, lengths)
    tests = compare_to_background(allc, ts, bg).merge(st_tests.drop(columns=["n_targets_bfd", "reason"]), on="tf_name")
    tests = tests.merge(tinfo[["tf_name", "tf_family", "tf_jgi_id"]], on="tf_name", how="left")
    stats = caveat_stats(ts, tests, allc, capped)
    lines = caveat_lines(stats)
    cl = mapped.merge(allc, on="bfd_protein_id", how="left")
    _write_commented(tables_tf / "T21_target_conservation.tsv", lines + [
        f"min_cov={MIN_COV} evalue={EVALUE} cov_mode=2 (query); n_orders/n_phyla/n_genera = distinct ORDER/PHYLUM/GENUS of genus representatives with a hit; "
        "n_genera also equals the number of representative genomes with a hit (one genome per genus)",
        f"background = all {len(bg)} FBF4BF29 proteins, search max-seqs {MAX_SEQS}"],
        cl[["tf_name", "tf_family", "jgi_protein_id", "bfd_protein_id", "n_phyla", "n_orders", "n_genera"]])
    _write_commented(tables_tf / "T21_target_background_tests.tsv", lines + [
        "p, p_adj: two-sided Mann-Whitney U of the TF's mapped targets against all FBF4BF29 proteins (n_orders), BH over tested TFs. "
        "effect: rank-biserial correlation, positive when targets are more widely conserved.",
        "effect_strat, p_strat, p_adj_strat: length-stratified (quintiles of the background protein length): effect = mean of the per-stratum rank-biserial "
        "effects over strata with >= 5 targets, weighted by the number of targets; p_strat = van Elteren stratified rank-sum test; BH over tested TFs. "
        "expected_median_n_orders: length-matched expected median (background reweighted to the TF's length profile). Per-stratum values: T21_target_strata.tsv."],
        tests)
    _write_commented(tables_tf / "T21_target_strata.tsv", lines[:1], strat)
    ecdf = _ecdf_table(allc, cl)
    ecdf.to_csv(tables_tf / "T21_target_ecdf.tsv", sep="\t", index=False)
    notes = notes_table(names, targets, ss, stats)
    notes.to_csv(tables_tf / "T21_target_notes.tsv", sep="\t", index=False)
    # ---- profile
    pf, sp = load_profile_sources(runs_tables, mapped["bfd_protein_id"].unique(), scratch)
    profile_targets(ts, pf, sp).to_csv(tables_tf / "T21_target_profile.tsv", sep="\t", index=False)
    print(f"tables in {tables_tf}", flush=True)


def _write_commented(path, comments, df):
    tmp = Path(str(path) + ".tmp")
    with open(tmp, "w") as f:
        for c in comments:
            f.write(f"# {c}\n")
        df.to_csv(f, sep="\t", index=False)
    tmp.replace(path)


def _capped_queries(stage: Path):
    """Queries whose hit table reached max-seqs rows (a lower bound: MMseqs2 caps before the e-value and coverage filters)."""
    gz = stage / "conservation_hits.tsv.gz"
    if not gz.exists():
        return None
    q = pd.read_csv(gz, sep="\t", header=None, usecols=[0], dtype=str)[0]
    c = q.value_counts()
    return int((c >= MAX_SEQS).sum()), int(len(c))


def _ecdf_table(allc, cl):
    x = np.arange(0, int(max(allc["n_orders"].max(), cl["n_orders"].max())) + 1)
    pooled = cl.drop_duplicates("bfd_protein_id")["n_orders"].to_numpy()
    noz = cl[cl["tf_name"] != "Znf1"].drop_duplicates("bfd_protein_id")["n_orders"].to_numpy()

    def ecdf(v):
        return np.searchsorted(np.sort(v), x, side="right") / max(len(v), 1)
    return pd.DataFrame({"n_orders": x, "ecdf_targets": ecdf(pooled), "ecdf_targets_no_Znf1": ecdf(noz),
                         "ecdf_background": ecdf(allc["n_orders"].to_numpy())})


def caveat_stats(ts, tests, allc, capped):
    m = ts[ts["bfd_protein_id"] != ""].drop_duplicates(["tf_name", "bfd_protein_id"])
    tested = tests[tests["reason"] == ""]
    per_prot = m[m["tf_name"].isin(tested["tf_name"])].groupby("bfd_protein_id").size()
    zn = m[m["tf_name"] == "Znf1"]["bfd_protein_id"].nunique()
    return {"mean_tfs_per_protein": float(per_prot.mean()), "n_unique_proteins": int(len(per_prot)),
            "znf1_share_rows": float((m["tf_name"] == "Znf1").sum() / max(len(m), 1)),
            "znf1_share_unique": float(zn / max(m["bfd_protein_id"].nunique(), 1)), "n_tested": int(len(tested)),
            "n_sig": int((tested["p_adj"] < 0.05).sum()), "n_sig_strat": int((tested["p_adj_strat"] < 0.05).sum()),
            "capped": capped, "frac_zero_hit_bg": float((allc["n_orders"] == 0).mean())}


def caveat_lines(s):
    cap = "not computed" if s["capped"] is None else f"{s['capped'][0]} of {s['capped'][1]} queries ({s['capped'][0] / max(s['capped'][1], 1):.0%}) reached it (lower bound)"
    return [
        "HYPOTHESIS ONLY. n_orders = number of orders among the genus representatives with a homology hit (e <= 1e-5, query coverage >= 0.5). "
        "It measures breadth of homology (domain and paralog hits), not ortholog conservation.",
        f"The {s['n_tested']} tests are dependent: the mapped targets overlap across TFs (mean {s['mean_tfs_per_protein']:.1f} tested TFs per target protein, "
        f"{s['n_unique_proteins']} unique proteins). {s['n_sig']} of {s['n_tested']} significant is not that many independent findings.",
        f"Znf1 holds {s['znf1_share_rows']:.0%} of the TF-target rows and {s['znf1_share_unique']:.0%} of the unique mapped target proteins. "
        "The background (all FBF4BF29 proteins) contains the targets.",
        f"Search max-seqs is {MAX_SEQS} per query; {cap}. Capped queries can have too few hits and so too low n_orders.",
        "Unmapped targets (different gene model or assembly) are not evidence of lineage-specific genes."]


def notes_table(names, targets, ss, stats):
    miss = names[(names["success_frip5_any"] == 1) & ~names["jgi_id"].isin(set(targets["tf_jgi_id"]))]
    rows = [("missing_tf", f"{r['name']} (jgi_id {r['jgi_id']}) is a successful TF of the paper with no sheet in MOESM4: it has no target set here.") for _, r in miss.iterrows()]
    d = ss[ss["n_mapped"] != ss["n_bfd_distinct"]]
    for _, r in d.iterrows():
        rows.append(("tie_collapse", f"{r['tf_name']}: {int(r['n_mapped'])} JGI targets map to {int(r['n_bfd_distinct'])} distinct BFD proteins "
                                     "(several JGI models share one BFD partner), so n_mapped differs from n_targets_bfd in the tests table."))
    rows += [("header_format", "T21_target_conservation.tsv, T21_target_background_tests.tsv and T21_target_strata.tsv start with lines that begin with '#' "
                               "(thresholds and caveats); read them with comment='#'."),
             ("header_format", "jgi_to_bfd_*.tsv (staging dir) start with one '#' line: min_pident, min_cov, cov_mode, tie rule."),
             ("map_sets", "Main map: best models (filtered GFF3, 10,700) on both sides, 90% identity, 0.8 coverage. Sensitivities: 70/0.7 and 90/0.5. "
                          "Secondary: reverse search against all 258,092 models (all_models_reverse); its rate is not one-to-one."),
             ("stratified", "Length strata are quintiles of the background protein length; see the comment lines of T21_target_background_tests.tsv.")]
    return pd.DataFrame(rows, columns=["topic", "note"])


def main(argv):
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("step", choices=["prep", "map", "repdb", "conserve", "tables", "figures"])
    ap.add_argument("--stage", default=str(RUNS / "analysis/tf_targets"))
    ap.add_argument("--scratch", default=None)
    ap.add_argument("--summary-dir", default=str(RUNS / "analysis/summary_tables"))
    ap.add_argument("--tables", default=str(RUNS / "tables"), help="BFD parquet tables (read only)")
    a = ap.parse_args(argv)
    stage = Path(a.stage)
    scratch = Path(a.scratch) if a.scratch else None
    if a.step == "prep":
        step_prep(stage)
    elif a.step == "map":
        step_map(stage, scratch)
    elif a.step == "repdb":
        step_repdb(stage, a.summary_dir)
    elif a.step == "conserve":
        step_conserve(stage, scratch, a.summary_dir)
    elif a.step == "tables":
        step_tables(stage, stage / "tables_tf", scratch, Path(a.tables))
    else:
        page_tf_targets(stage, stage / "tables_tf")


if __name__ == "__main__":
    main(sys.argv[1:])
