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
                      split_memory_limit="30G", max_seqs=3000, sensitivity=4.0):
    """MMseqs2 search of query_fasta against the genus-representative database. Query coverage is
    --cov-mode 2. max_seqs caps the hits per query; sensitivity is below the default 5.7 to fit the time limit."""
    exe = jm.need_mmseqs()
    tmp = Path(tmp)
    tmp.mkdir(parents=True, exist_ok=True)
    cmd = [exe, "easy-search", str(query_fasta), str(db), str(hits_tsv), str(tmp / "mm_tmp"),
           "-e", str(evalue), "-c", str(min_cov), "--cov-mode", "2", "-s", str(sensitivity),
           "--max-seqs", str(max_seqs), "--format-output", "query,target,evalue,qcov,bits",
           "--threads", str(threads or jm.threads_default()), "--split-memory-limit", str(split_memory_limit), "-v", "1"]
    subprocess.run(cmd, check=True)


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
               "background_median_n_orders": float(np.median(bg)), "effect": np.nan, "p": np.nan, "reason": ""}
        if len(vals) < MIN_TARGETS:
            row["reason"] = f"fewer than {MIN_TARGETS} mapped targets"
        else:
            u, p = st.mannwhitneyu(vals, bg, alternative="two-sided")
            row["effect"] = float(2 * u / (len(vals) * len(bg)) - 1)
            row["p"] = float(p)
        rows.append(row)
    out = pd.DataFrame(rows, columns=["tf_name", "n_targets_bfd", "median_n_orders", "background_median_n_orders",
                                      "effect", "p", "reason"])
    out["p_adj"] = bh_adjust(out["p"].to_numpy())
    return out[["tf_name", "n_targets_bfd", "median_n_orders", "background_median_n_orders", "effect", "p", "p_adj", "reason"]]


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
def family_colors(families):
    import matplotlib.pyplot as plt
    pal = list(plt.get_cmap("tab20").colors) + list(plt.get_cmap("tab20b").colors)
    return {f: pal[i % len(pal)] for i, f in enumerate(sorted(set(families)))}


def _fig_conservation(cons, tests):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    ok = tests[tests["reason"].fillna("") == ""].sort_values("median_n_orders")
    if ok.empty:
        return None
    fam = cons.drop_duplicates("tf_name").set_index("tf_name")["tf_family"]
    cm = family_colors(fam[ok["tf_name"]])
    fig, ax = plt.subplots(figsize=(10.0, max(4.0, 0.2 * len(ok) + 1.6)))
    data = [cons.loc[cons["tf_name"] == t, "n_orders"].to_numpy() for t in ok["tf_name"]]
    bp = ax.boxplot(data, orientation="horizontal", patch_artist=True, widths=0.7, showfliers=False)
    for patch, t in zip(bp["boxes"], ok["tf_name"]):
        patch.set_facecolor(cm[fam[t]])
        patch.set_alpha(0.9)
    for med in bp["medians"]:
        med.set_color("black")
    bgm = float(ok["background_median_n_orders"].iloc[0])
    ax.axvline(bgm, color="#c0392b", ls="--", lw=1.2)
    ax.set_yticks(range(1, len(ok) + 1))
    ax.set_yticklabels([f"{t} (n={n})" for t, n in zip(ok["tf_name"], ok["n_targets_bfd"])], fontsize=6)
    ax.set_xlabel("Orders with a hit (genus representatives; e <= 1e-5, query coverage >= 0.5)")
    ax.set_title(f"Conservation of DAP-seq targets (mapped, n >= 20); dashed line = median of all proteins ({bgm:g})", fontsize=8)
    handles = [Patch(color=c, label=f) for f, c in cm.items()]
    ax.legend(handles=handles, fontsize=5.5, title="paper TF family", title_fontsize=6.5, loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False)
    fig.tight_layout()
    return fig


def _fig_sizes(sizes):
    import matplotlib.pyplot as plt
    s = sizes.sort_values("n_targets")
    fig, ax = plt.subplots(figsize=(7.5, max(4.0, 0.15 * len(s) + 1.6)))
    y = np.arange(len(s))
    ax.barh(y, s["n_mapped"], color="#2a6f97", label="mapped to BFD (RBH)")
    ax.barh(y, s["n_unmapped"], left=s["n_mapped"], color="#c9c9c9", label="no BFD partner")
    ax.set_yticks(y)
    ax.set_ylim(-0.7, len(s) - 0.3)
    ax.set_yticklabels(s["tf_name"], fontsize=5)
    ax.set_xlabel("DAP-seq targets (JGI protein IDs)")
    ax.set_title("Targets with and without a BFD partner, per TF", fontsize=9)
    ax.legend(fontsize=7, frameon=False, loc="lower right")
    fig.tight_layout()
    return fig


def page_tf_targets(out_dir, tables_dir=None, o=None):
    """Draw T21_target_conservation and T21_target_set_sizes from the T21 tables. tables_dir default:
    <out_dir>/tables_tf. Figures go to <out_dir>/figures/tf (Out level 'tf'). Returns the list of Made records."""
    from figs_common import safe
    from make_figures import Out
    out_dir = Path(out_dir)
    td = Path(tables_dir) if tables_dir else out_dir / "tables_tf"
    o = o or Out(out_dir, "tf")
    cons = pd.read_csv(td / "T21_target_conservation.tsv", sep="\t", comment="#")
    tests = pd.read_csv(td / "T21_target_background_tests.tsv", sep="\t", keep_default_na=False, na_values=[""])
    tests["reason"] = tests["reason"].fillna("")
    sizes = pd.read_csv(td / "T21_tf_set_sizes.tsv", sep="\t")
    ntest = int((tests["reason"] == "").sum())
    o.info = {"n_genomes": None, "n_groups": ntest, "pooled_share": 0.0}
    o.meta["T21_target_conservation"] = {"signal": ntest, "signal_kind": "n_tfs_tested"}
    safe(o, "T21_target_conservation", "tf", lambda: _fig_conservation(cons, tests))
    safe(o, "T21_target_set_sizes", "tf", lambda: _fig_sizes(sizes))
    return o.made


# ---------------------------------------------------------------- pipeline steps
def _gz_write(df, path):
    with gzip.open(path, "wt") as f:
        df.to_csv(f, sep="\t", index=False)


def step_prep(stage: Path):
    """Write the FBF4BF29 proteome and the JGI protein FASTA (all models) with plain ID headers."""
    stage.mkdir(parents=True, exist_ok=True)
    bfd = stage / "FBF4BF29.faa"
    if not bfd.exists():
        src = proteome_paths(RUNS / "seqdb/genomes.tsv")[BFD_GENOME]
        tmp = bfd.with_suffix(".tmp")
        n = 0
        with jm._open_text(src) as f, open(tmp, "w") as o:
            for line in f:
                if line.startswith(">"):
                    o.write(f">{line[1:].split()[0]}\n")
                    n += 1
                else:
                    o.write(line)
        tmp.rename(bfd)
        print(f"wrote {bfd} ({n} proteins) from {src}", flush=True)
    jf = stage / "Rhimi59_2_all_models.faa"
    if not jf.exists():
        prot = jm.read_jgi_proteins(JGI_FASTA)
        tmp = jf.with_suffix(".tmp")
        with open(tmp, "w") as o:
            for pid, (_, seq) in prot.items():
                o.write(f">{pid}\n{seq}\n")
        tmp.rename(jf)
        print(f"wrote {jf} ({len(prot)} proteins)", flush=True)


def step_map(stage: Path, scratch: Path):
    """RBH map at 90/0.8 (default) and 70/0.7 (sensitivity); then the rate table."""
    for tag, pid, cov in (("p90_c0.8", 90.0, 0.8), ("p70_c0.7", 70.0, 0.7)):
        out = stage / f"jgi_to_bfd_{tag}.tsv"
        if out.exists():
            continue
        tmp = out.with_suffix(".tmp")
        jm.map_jgi_to_bfd(stage / "Rhimi59_2_all_models.faa", stage / "FBF4BF29.faa", tmp, scratch / f"map_{tag}", pid, cov)
        tmp.rename(out)
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
    catalog = jm.read_jgi_catalog(JGI_GFF)
    best = set(catalog["protein_id"])
    # map rates (default and sensitivity thresholds)
    rates, maps = [], {}
    uniq_targets = targets["jgi_protein_id"].unique()
    tested_ids = tested["jgi_protein_id"].unique()
    for tag in ("p90_c0.8", "p70_c0.7"):
        m = jm.read_map(stage / f"jgi_to_bfd_{tag}.tsv")
        maps[tag] = m
        for what, ids in (("tested_tfs", tested_ids), ("targets", uniq_targets)):
            r = jm.mapping_rate(m, ids, best_ids=best)
            rates.append({"thresholds": tag, "set": what, **r})
        r = jm.mapping_rate(m, list(m["jgi_protein_id"]), best_ids=best)
        rates.append({"thresholds": tag, "set": "all_jgi_models", **r})
    pd.DataFrame(rates).to_csv(tables_tf / "T21_jgi_map_rate.tsv", sep="\t", index=False)
    ts = target_sets(targets, maps["p90_c0.8"])
    ts = ts.merge(targets[["tf_name", "tf_family"]].drop_duplicates("tf_name"), on="tf_name", how="left")
    ts["in_best_models"] = ts["jgi_protein_id"].isin(best)
    _gz_write(ts, tables_tf / "T21_target_sets.tsv.gz")
    ss = set_sizes(ts)
    tinfo = targets.drop_duplicates("tf_name")[["tf_name", "tf_family", "tf_jgi_id"]].merge(
        names[["jgi_id", "success_frip5_any"]], left_on="tf_jgi_id", right_on="jgi_id", how="left").drop(columns=["jgi_id"])
    ss = ss.merge(tinfo, on="tf_name", how="left")
    key = tinfo
    ss.to_csv(tables_tf / "T21_tf_set_sizes.tsv", sep="\t", index=False)
    # conservation
    allc = pd.read_csv(stage / "T21_conservation_all_bfd.tsv", sep="\t", comment="#")
    bg = jm.fasta_ids(stage / "FBF4BF29.faa")
    mapped = ts[ts["bfd_protein_id"] != ""]
    cl = mapped.merge(allc, on="bfd_protein_id", how="left")
    with open(tables_tf / "T21_target_conservation.tsv", "w") as f:
        f.write(f"# min_cov={MIN_COV} evalue={EVALUE} cov_mode=2 n_orders/n_phyla/n_genera of genus representatives; background = all {len(bg)} FBF4BF29 proteins\n")
        cl[["tf_name", "tf_family", "jgi_protein_id", "bfd_protein_id", "n_phyla", "n_orders", "n_genera"]].to_csv(f, sep="\t", index=False)
    tests = compare_to_background(allc, ts, bg)
    tests.to_csv(tables_tf / "T21_target_background_tests.tsv", sep="\t", index=False)
    # profile
    pf, sp = load_profile_sources(runs_tables, mapped["bfd_protein_id"].unique(), scratch)
    profile_targets(ts, pf, sp).to_csv(tables_tf / "T21_target_profile.tsv", sep="\t", index=False)
    print(f"tables in {tables_tf}", flush=True)


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
