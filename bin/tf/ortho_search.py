"""Ortholog search for WC-2, Zap1-like and bZip5 among the genus-representative proteomes.

Steps (each skips when its output exists; run from the repo root):
    python bin/tf/ortho_search.py seeds       seed table, seed FASTA per target, best BFD match of every seed
    python bin/tf/ortho_search.py candidates  proteins of the genus representatives with the target Pfam, with architecture
    python bin/tf/ortho_search.py search      MMseqs2 seeds vs candidates, then the architecture rule

Targets: WC-2 (GATA and exactly one merged PAS region), WC-1 (GATA and two or more merged PAS regions), Zap1-like
(C2H2 >= 4), bZip5 (bZIP_1). WC-1 and WC-2 share one candidate set (all GATA proteins) and one search; each candidate is
assigned to its best seed (highest bit score) over the WC-1 and WC-2 seeds. A candidate whose architecture does not agree
with the rule of its best seed is not counted and is reported in the disagreement column.
Coverage: MMseqs2 --cov-mode 2 (>= min_cov of the SEED sequence inside the alignment). Mode 0 (both sequences) kept only
proteins of the seed length (ZAP1 880 aa: 53 genera). Mode 5 (shorter sequence) was tried and rejected: it lets any
candidate that contains a short seed-like stretch pass, so 10610 of 14301 GATA proteins and 29799 of 45364 C2H2 proteins hit.
Mode 2 needs half of the seed to align, so a short or truncated candidate can pass and a long candidate is not penalised.
It still admits paralogs (C2H2 repeats); the architecture rule and the tree limit that.
Order of the rules: candidates carry the target Pfam only. The search (identity and coverage thresholds) runs on all of
them. The architecture rule (for example GATA plus PAS) is applied to the hits AFTER the search. A genus
representative without the target domain is kept in the candidate table as an absent row, not dropped.
Similarity and tree position are evidence for orthology, not proof.
"""
import argparse
import gzip
import hashlib
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "bin" / "summary"))
import jgi_map as jm  # noqa: E402
import seeds as sd  # noqa: E402

RUNS = Path("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs")
JGI_FASTA = (Path("/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/resources/jgi_atcc_11559/Rhimi59_2/Mycocosm/Annotation")
             / "All_models,_Filtered_and_Not/Proteins/Rhimi59_2_all_proteins_20220409.aa.fasta.gz")
MIN_PIDENT = 30.0
MIN_COV = 0.5
COV_MODE = 2
SLUG = {"WC-2": "wc2", "WC-1": "wc1", "Zap1-like": "zap1", "bZip5": "bzip5"}
GROUP = {"WC-2": "wc", "WC-1": "wc", "Zap1-like": "zap1", "bZip5": "bzip5"}      # targets that share candidates and search
PAS_ACCS = sd.PAS_ACCS
pas_regions = sd.pas_regions
DISTANT_NOTE = "no BFD FBF4BF29 match at 30%/0.5 (distant seed; not an error)"
ANCHORS = {"FBF4BF29": "Rhizopus microsporus ATCC 11559", "F750FE6D": "Mucor lusitanicus", "F43EADC6": "Phycomyces blakesleeanus"}
BFD_GENOME = "FBF4BF29"
HIT_COLS = jm.HIT_COLS


# ---------------------------------------------------------------- rules
def _cmp(k, op, n):
    return k == n if op == "==" else k >= n


def architecture_ok(arch: str, required, n_pas=None) -> bool:
    """True when the architecture string (accessions joined by '>', no versions) meets every entry of `required`.
    Entry syntax: see seeds.parse_rule. n_pas is the number of merged PAS regions of the protein (needed by
    PAS_regions rules; a PAS_regions rule fails when it is None). An empty list is always true."""
    doms = [d for d in str(arch).split(">") if d] if isinstance(arch, str) else []
    for entry in required:
        r = sd.parse_rule(entry)
        if r[0] == "any":
            ok = bool(r[1] & set(doms))
        elif r[0] == "count":
            ok = _cmp(sum(1 for d in doms if d in r[1]), r[2], r[3])
        else:
            ok = n_pas is not None and _cmp(int(n_pas), r[1], r[2])
        if not ok:
            return False
    return True


def architecture_names(arch, names):
    """Accession architecture -> readable text with Pfam names (names: accession -> name)."""
    return ">".join(names.get(a, a) for a in str(arch).split(">") if a)


# ---------------------------------------------------------------- candidates
def select_candidates(genus_reps: pd.DataFrame, pfam_hits: pd.DataFrame, target: str, seeds: pd.DataFrame) -> pd.DataFrame:
    """All proteins of the genus representatives that carry the target Pfam (seeds['target_pfam'] of `target`).
    pfam_hits: every Pfam domain row of those proteins (protein_id, pfam_acc, ali_from, ali_to). Columns of the result:
    species_prefix, protein_id, architecture, n_target_domains, n_pas_regions (merged PAS regions), gata_rel_pos (start of
    the first target-Pfam domain divided by the end of the last domain), present. A representative with no such protein
    gets one row with protein_id NaN, empty architecture and present False (reported as absent)."""
    tpf = seeds.loc[seeds["target"] == target, "target_pfam"].iloc[0]
    ph = pfam_hits.copy()
    ph["acc"] = ph["pfam_acc"].map(sd.strip_version)
    if "ali_to" not in ph.columns:
        ph["ali_to"] = ph["ali_from"]
    ph = ph.sort_values(["protein_id", "ali_from"], kind="stable")
    arch = ph.groupby("protein_id")["acc"].agg(">".join)
    ntar = ph[ph["acc"] == tpf].groupby("protein_id").size()
    first = ph[ph["acc"] == tpf].groupby("protein_id")["ali_from"].min()
    last = ph.groupby("protein_id")["ali_to"].max()
    pas = ph[ph["acc"].isin(PAS_ACCS)]
    npas = {pid: len(pas_regions(g)) for pid, g in pas.groupby("protein_id")}
    reps = set(genus_reps["species_prefix"])
    rows = [{"species_prefix": str(p).split("_")[0], "protein_id": p, "architecture": arch[p],
             "n_target_domains": int(n), "n_pas_regions": int(npas.get(p, 0)),
             "gata_rel_pos": float(first[p] / last[p]) if last[p] else np.nan, "present": True}
            for p, n in ntar.items() if str(p).split("_")[0] in reps]
    cols = ["species_prefix", "protein_id", "architecture", "n_target_domains", "n_pas_regions", "gata_rel_pos", "present"]
    df = pd.DataFrame(rows, columns=cols)
    have = set(df["species_prefix"])
    absent = pd.DataFrame([{"species_prefix": g, "protein_id": np.nan, "architecture": "", "n_target_domains": 0,
                            "n_pas_regions": 0, "gata_rel_pos": np.nan, "present": False}
                           for g in genus_reps["species_prefix"] if g not in have], columns=cols)
    out = pd.concat([df, absent], ignore_index=True)
    return out.sort_values(["species_prefix", "protein_id"], na_position="last", kind="stable").reset_index(drop=True)


# ---------------------------------------------------------------- search
def write_commented(df: pd.DataFrame, path, params: dict):
    """Write a TSV whose first line is '# key=value ...' with the thresholds. Read it with comment='#'."""
    with open(path, "w") as f:
        f.write("# " + " ".join(f"{k}={v}" for k, v in params.items()) + "\n")
        df.to_csv(f, sep="\t", index=False)


def search(seed_fasta: Path, cand_fasta: Path, out: Path, scratch: Path, min_pident=MIN_PIDENT, min_cov=MIN_COV,
           threads=None, cov_mode=COV_MODE) -> pd.DataFrame:
    """MMseqs2 easy-search of the seeds against the candidates. Thresholds: identity min_pident (percent), coverage
    min_cov of the seed (cov mode 2). The thresholds are written as a '# ' line at the top of `out` and kept in
    df.attrs. Without mmseqs on PATH the step is skipped: an empty table, attrs['skipped'] says why."""
    params = {"min_pident": float(min_pident), "min_cov": float(min_cov), "cov_mode": cov_mode, "max_seqs": 100000}
    if shutil.which("mmseqs") is None:
        msg = "mmseqs is not on PATH: search skipped (module load MMseqs2/17-b804f)"
        print(msg, flush=True)
        df = pd.DataFrame(columns=HIT_COLS)
        df.attrs.update(params, skipped=msg)
        return df
    Path(scratch).mkdir(parents=True, exist_ok=True)
    raw = Path(str(out) + ".raw")
    df = jm.easy_search(Path(seed_fasta), Path(cand_fasta), raw, Path(scratch), min_pident, min_cov, threads=threads,
                        max_seqs=params["max_seqs"], cov_mode=cov_mode)
    raw.unlink(missing_ok=True)
    write_commented(df, out, params)
    df.attrs.update(params)
    return df


def pass_hits(hits: pd.DataFrame, cand: pd.DataFrame, required) -> pd.DataFrame:
    """Candidates with a hit that passed the search thresholds AND the architecture rule. One row per candidate protein,
    with its best seed (highest bit score). Columns: species_prefix, protein_id, architecture, seed, pident, qcov,
    tcov, evalue, bits."""
    cols = ["species_prefix", "protein_id", "architecture", "seed", "pident", "qcov", "tcov", "evalue", "bits"]
    if hits is None or not len(hits):
        return pd.DataFrame(columns=cols)
    best = hits.sort_values("bits", ascending=False, kind="stable").drop_duplicates("target")
    ccols = [c for c in ("species_prefix", "protein_id", "architecture", "n_pas_regions") if c in cand.columns]
    m = best.merge(cand[cand["present"]][ccols], left_on="target", right_on="protein_id", how="inner")
    if "n_pas_regions" not in m.columns:
        m["n_pas_regions"] = None
    m = m[[architecture_ok(a, required, n) for a, n in zip(m["architecture"], m["n_pas_regions"])]]
    m = m.rename(columns={"query": "seed"})
    return m[cols].sort_values(["species_prefix", "protein_id"]).reset_index(drop=True)


def assign_best_seed(hits: pd.DataFrame, cand: pd.DataFrame, seed_target: dict, rules: dict) -> pd.DataFrame:
    """One row per candidate with a hit: its best seed over ALL seeds of the group (highest bit score), the target of that
    seed (seed_target), and whether the architecture rule of that target holds (rule_ok). arch_type names the target
    whose rule the protein meets (first match in `rules` order; empty if none). disagree is True when rule_ok is False:
    these rows are reported, not dropped. rules: target -> list of required entries."""
    cols = ["species_prefix", "protein_id", "architecture", "n_pas_regions", "seed", "seed_target", "pident", "qcov", "tcov",
            "evalue", "bits", "rule_ok", "arch_type", "disagree"]
    if hits is None or not len(hits):
        return pd.DataFrame(columns=cols)
    best = hits.sort_values("bits", ascending=False, kind="stable").drop_duplicates("target")
    ccols = [c for c in ("species_prefix", "protein_id", "architecture", "n_pas_regions") if c in cand.columns]
    m = best.merge(cand[cand["present"]][ccols], left_on="target", right_on="protein_id", how="inner")
    if "n_pas_regions" not in m.columns:
        m["n_pas_regions"] = None
    m = m.rename(columns={"query": "seed"})
    m["seed_target"] = m["seed"].map(seed_target)
    ok, typ = [], []
    for a, n, t in zip(m["architecture"], m["n_pas_regions"], m["seed_target"]):
        ok.append(architecture_ok(a, rules.get(t, []), n))
        typ.append(next((tt for tt, rq in rules.items() if architecture_ok(a, rq, n)), ""))
    m["rule_ok"], m["arch_type"] = ok, typ
    m["disagree"] = ~m["rule_ok"]
    return m[cols].sort_values(["species_prefix", "protein_id"]).reset_index(drop=True)


def pass_for_target(assigned: pd.DataFrame, target: str) -> pd.DataFrame:
    """Proteins whose best seed is a seed of `target` and whose architecture agrees with the rule of `target`."""
    return assigned[(assigned["seed_target"] == target) & assigned["rule_ok"].astype(bool)].reset_index(drop=True)


def count_disagreements(assigned: pd.DataFrame, target: str) -> int:
    """Proteins whose best seed is a seed of `target` but whose architecture fails the rule of `target`."""
    return int(((assigned["seed_target"] == target) & ~assigned["rule_ok"].astype(bool)).sum())


# ---------------------------------------------------------------- re-run logic
STAMP_NAME = ".config_stamp"
OUTPUT_GLOBS = ["seeds.tsv", "seeds_*.faa", "cand_*", "hits_*", "pass_*", "assigned_*", "tree_*", "counts.tsv"]


def config_hash(seed_config: Path, code_files, **thresholds) -> str:
    """Hash of the seed config, the thresholds and the code files. A change means the outputs are stale."""
    h = hashlib.sha256()
    for f in [Path(seed_config)] + [Path(c) for c in code_files]:
        h.update(f.read_bytes() if f.exists() else b"<absent>")
    th = dict(min_pident=MIN_PIDENT, min_cov=MIN_COV, cov_mode=COV_MODE, pas=",".join(PAS_ACCS))
    th.update(thresholds)
    h.update(repr(sorted(th.items())).encode())
    return h.hexdigest()


def invalidate_if_changed(stage: Path, root: Path, code_files=None, **thresholds) -> bool:
    """Compare the stored stamp with the current hash. On a change (or no stamp) remove the derived outputs of the stage
    (seeds, candidates, hits, passes, trees, counts; not the Pfam cache or the UniProt seed cache), write the new stamp, and
    return True. Return False when nothing changed."""
    code_files = [Path(__file__), HERE / "seeds.py", HERE / "ortho_tree.py"] if code_files is None else code_files
    cur = config_hash(Path(root) / "config" / "tf_seeds.tsv", code_files, **thresholds)
    sf = Path(stage) / STAMP_NAME
    if sf.exists() and sf.read_text().strip() == cur:
        return False
    for g in OUTPUT_GLOBS:
        for f in Path(stage).glob(g):
            if f.is_file() or f.is_symlink():
                f.unlink()
    sf.write_text(cur + "\n")
    return True


def link_reuse(link: Path, target: Path):
    """Symlink `link` to `target`; replaces a dangling or stale symlink, leaves a good one alone."""
    link = Path(link)
    if link.is_symlink() and not link.exists():
        link.unlink()
    if not (link.is_symlink() or link.exists()):
        link.symlink_to(target)


# ---------------------------------------------------------------- data access
def read_pfam_rows(tables_dir, reps, accs):
    """Every Pfam domain row of the proteins that are in the representative genomes and carry one of `accs`.
    Reads tables/pfam.parquet (read only). The protein list is the genome prefix list of the representatives."""
    import duckdb
    con = duckdb.connect()
    plist = ",".join("'" + g + "'" for g in sorted(reps))
    alist = ",".join("'" + a + "'" for a in accs)
    pq = f"'{Path(tables_dir) / 'pfam.parquet'}'"
    return con.execute(f"""
        with t as (select distinct protein_id from {pq}
                   where regexp_extract(pfam_acc, '^(PF[0-9]+)', 1) in ({alist})
                     and split_part(protein_id, '_', 1) in ({plist}))
        select p.protein_id, p.pfam_id, p.pfam_acc, p.ali_from, p.ali_to, p.domain_i_evalue
        from {pq} p join t using (protein_id)""").df()


def read_pfam_for_proteins(tables_dir, ids):
    import duckdb
    con = duckdb.connect()
    con.register("ids", pd.DataFrame({"protein_id": list(ids)}))
    pq = f"'{Path(tables_dir) / 'pfam.parquet'}'"
    return con.execute(f"select p.protein_id, p.pfam_id, p.pfam_acc, p.ali_from from {pq} p join ids using (protein_id)").df()


def extract_from_db(db: Path, ids, out_fa, scratch: Path):
    """Write the sequences of `ids` from the genus-representative MMseqs2 database (createsubdb, convert2fasta).
    Returns the number written. The header is cut to the first word."""
    import duckdb
    exe = jm.need_mmseqs()
    scratch = Path(scratch)
    scratch.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect()
    con.register("ids", pd.DataFrame({"name": sorted(set(ids))}))
    keys = con.execute(f"""select cast(k as bigint) from read_csv('{db}.lookup', delim='\t', header=false,
                               columns={{'k': 'VARCHAR', 'name': 'VARCHAR', 'f': 'VARCHAR'}}) l join ids using (name)""").fetchall()
    kf = scratch / "keys.txt"
    kf.write_text("".join(f"{k[0]}\n" for k in keys))
    sub = scratch / "sub"
    subprocess.run([exe, "createsubdb", str(kf), str(db), str(sub), "-v", "1"], check=True)
    subprocess.run([exe, "createsubdb", str(kf), str(db) + "_h", str(sub) + "_h", "-v", "1"], check=True)
    subprocess.run([exe, "convert2fasta", str(sub), str(scratch / "sub.fa"), "-v", "1"], check=True)
    n = 0
    with open(scratch / "sub.fa") as f, gzip.open(out_fa, "wt") if str(out_fa).endswith(".gz") else open(out_fa, "w") as o:
        for line in f:
            if line.startswith(">"):
                o.write(">" + line[1:].split()[0] + "\n")
                n += 1
            else:
                o.write(line)
    for p in scratch.glob("sub*"):
        p.unlink()
    return n


def filter_fasta(src, ids, dst):
    ids = set(ids)
    keep, n = False, 0
    with jm._open_text(src) as f, (gzip.open(dst, "wt") if str(dst).endswith(".gz") else open(dst, "w")) as o:
        for line in f:
            if line.startswith(">"):
                keep = line[1:].split()[0] in ids
                n += keep
            if keep:
                o.write(line)
    return n


# ---------------------------------------------------------------- pipeline steps
def load_cfg(root):
    return sd.read_seed_config(Path(root) / "config" / "tf_seeds.tsv")


def arch_map(pf):
    ph = pf.assign(acc=pf["pfam_acc"].map(sd.strip_version)).sort_values(["protein_id", "ali_from"], kind="stable")
    return ph.groupby("protein_id")["acc"].agg(">".join).to_dict()


def npas_map(pf):
    pas = pf[pf["pfam_acc"].map(sd.strip_version).isin(PAS_ACCS)]
    return {pid: len(pas_regions(g)) for pid, g in pas.groupby("protein_id")}


def ROOT_NAMES(root):
    return Path(root) / "config" / "paper_tf_names.tsv"


def jm_proteome_paths():
    from targets import proteome_paths
    return proteome_paths(RUNS / "seqdb/genomes.tsv")


def read_pfam_for_genomes(tables_dir, genomes):
    import duckdb
    con = duckdb.connect()
    gl = ",".join("'" + g + "'" for g in genomes)
    pq = f"'{Path(tables_dir) / 'pfam.parquet'}'"
    return con.execute(f"select protein_id, pfam_id, pfam_acc, ali_from, ali_to from {pq} where split_part(protein_id, '_', 1) in ({gl})").df()


def step_seeds(stage: Path, root: Path, tables_dir: Path, scratch: Path):
    out = stage / "seeds.tsv"
    if out.exists():
        print(f"seeds: {out} exists, skipped", flush=True)
        return pd.read_csv(out, sep="\t", keep_default_na=False)
    cfg = load_cfg(root)
    rows = []
    pap = sd.paper_sequences(cfg, JGI_FASTA)          # paper seeds: full-length JGI sequences
    for (_, c), (_, p) in zip(cfg[cfg["source"] == "paper"].iterrows(), pap.iterrows()):
        rows.append({"target": c["target"], "label": c["seed_name"], "source": "paper", "key": c["key"], "seed_id": p["accession"],
                     "accession": p["accession"], "sequence": p["sequence"], "status": p["status"], "genome": BFD_GENOME,
                     "phylum": "Mucoromycota"})
    uni = sd.fetch_seeds(cfg, stage / "seed_cache")   # UniProt seeds, cached
    for (_, c), (_, u) in zip(cfg[cfg["source"] == "uniprot"].iterrows(), uni.iterrows()):
        rows.append({"target": c["target"], "label": c["seed_name"], "source": "uniprot", "key": c["key"],
                     "seed_id": u["accession"] or f"uniprot_{c['key']}", "accession": u["accession"], "sequence": u["sequence"],
                     "status": u["status"], "genome": "", "phylum": sd.TAXON_PHYLUM.get(c["taxon_id"], "")})
    bfd_rows = cfg[cfg["source"] == "bfd"]              # BFD seeds selected by domain rule
    if len(bfd_rows):
        paths = jm_proteome_paths()
        pf = read_pfam_for_genomes(tables_dir, sorted({k.split("|")[0] for k in bfd_rows["key"]}))
        seqcache = {}
        for _, c in bfd_rows.iterrows():
            g, sel = c["key"].split("|")
            ids = sd.bfd_selector_proteins(g, sel, pf)
            if g not in seqcache:
                seqcache[g] = sd.read_fasta_dict(paths[g]) if g in paths and Path(paths[g]).exists() else {}
            if not ids:
                rows.append({"target": c["target"], "label": c["seed_name"], "source": "bfd", "key": c["key"], "seed_id": "",
                             "accession": "", "sequence": "", "status": f"missing: no protein of {g} matches {sel}",
                             "genome": g, "phylum": "Mucoromycota"})
            for pid in ids:
                rows.append({"target": c["target"], "label": f"{c['seed_name']} {pid}", "source": "bfd", "key": c["key"],
                             "seed_id": pid, "accession": pid, "sequence": seqcache[g].get(pid, ""),
                             "status": "ok" if seqcache[g].get(pid) else "missing: protein sequence not found", "genome": g,
                             "phylum": "Mucoromycota"})
    sdf = pd.DataFrame(rows)
    # best BFD (FBF4BF29) match of every seed with the thresholds of the ortholog search, and its architecture
    bfd_faa = stage / "FBF4BF29.faa"
    ok = sdf[sdf["status"] == "ok"]
    q = scratch / "all_seeds.faa"
    q.write_text("".join(f">{r.seed_id}\n{r.sequence}\n" for r in ok.itertuples()))
    hits = jm.easy_search(q, bfd_faa, scratch / "seeds_vs_bfd.tsv", scratch / "mm", MIN_PIDENT, MIN_COV, max_seqs=20, cov_mode=COV_MODE)
    best = hits.sort_values("bits", ascending=False).drop_duplicates("query").set_index("query") if len(hits) else pd.DataFrame()
    pfb = read_pfam_for_genomes(tables_dir, [BFD_GENOME])
    arch, npas = arch_map(pfb), npas_map(pfb)
    names = pd.read_parquet(RUNS / "analysis/summary_tables/pfam_names.parquet").set_index("pfam_acc")["pfam_id"].to_dict()
    mt, mp, mq, ma, man, mn = [], [], [], [], [], []
    for r in sdf.itertuples():
        if r.seed_id in best.index:
            b = best.loc[r.seed_id]
            a = arch.get(b["target"], "")
            mt.append(b["target"]); mp.append(float(b["pident"])); mq.append(float(b["qcov"])); ma.append(a)
            man.append(architecture_names(a, names)); mn.append(int(npas.get(b["target"], 0)))
        else:
            mt.append(""); mp.append(np.nan); mq.append(np.nan); ma.append(""); man.append(""); mn.append(0)
    sdf["bfd_match"], sdf["bfd_pident"], sdf["bfd_seed_cov"] = mt, mp, mq
    sdf["bfd_architecture"], sdf["bfd_architecture_names"], sdf["bfd_n_pas_regions"] = ma, man, mn
    # paper-named sanity check: resolve against BFD with the paper names (thresholds as the Task 5b JGI-to-BFD map)
    pap_seqs = {r.seed_id: r.sequence for r in sdf[(sdf.source == "paper") & (sdf.status == "ok")].itertuples()}
    res = sd.resolve_paper_seeds(pap_seqs, bfd_faa, ROOT_NAMES(root), arch, scratch / "resolve", min_pident=90.0, min_cov=0.8)
    rmap = res.set_index("paper_id")
    sdf["paper_name"] = [rmap.loc[i, "name"] if i in rmap.index else "" for i in sdf["seed_id"]]
    sdf["warning"] = ""
    for i, r in sdf.iterrows():
        w = ""
        if r["source"] == "paper" and r["seed_id"] in rmap.index:
            w = rmap.loc[r["seed_id"], "warning"]
            sdf.loc[i, "bfd_match_paper_rule"] = rmap.loc[r["seed_id"], "bfd_protein"]
        elif r["source"] == "uniprot" and r["status"] == "ok":
            tpf, req = sd.rules_for_target(cfg, r["target"])
            if not r["bfd_match"]:
                w = DISTANT_NOTE
            elif not architecture_ok(r["bfd_architecture"], ([tpf] if tpf else []) + req, r["bfd_n_pas_regions"]):
                w = (f"architecture check failed: best BFD match {r['bfd_match']} has '{r['bfd_architecture']}' "
                     f"({r['bfd_n_pas_regions']} PAS regions); rule for {r['target']}: {[tpf] + req}")
        sdf.loc[i, "warning"] = w
    sdf.to_csv(out, sep="\t", index=False)
    for grp in sorted(set(GROUP.values())):
        tg = [t for t in GROUP if GROUP[t] == grp]
        g = sdf[(sdf["target"].isin(tg)) & (sdf["status"] == "ok")]
        (stage / f"seeds_{grp}.faa").write_text("".join(f">{r.seed_id}\n{r.sequence}\n" for r in g.itertuples()))
    for t in SLUG:
        sub = sdf[sdf["target"] == t]
        print(f"seeds {t}: {int((sub.status == 'ok').sum())} ok, {int((sub.status != 'ok').sum())} missing; "
              + "; ".join(f"{r.label}: {r.status}" for r in sub.itertuples() if r.source == "uniprot"), flush=True)
    return sdf


def step_candidates(stage: Path, root: Path, tables_dir: Path, scratch: Path):
    cfg = load_cfg(root)
    reps = pd.read_csv(stage / "genus_reps.tsv", sep="\t")
    accs = sorted({sd.strip_version(sd.rules_for_target(cfg, t)[0]) for t in SLUG})
    pf_path = stage / "pfam_candidates.parquet"
    if pf_path.exists():
        pf = pd.read_parquet(pf_path)
    else:
        pf = read_pfam_rows(tables_dir, reps["species_prefix"], accs)
        pf.to_parquet(pf_path, index=False)
    print(f"pfam rows of candidate proteins: {len(pf)}", flush=True)
    groups = {}
    for t, g in GROUP.items():
        groups.setdefault(g, t)                      # one candidate set per group (targets of a group share the target Pfam)
    todo = [g for g in groups if not (stage / f"cand_{g}.faa.gz").exists()]
    tabs, union_ids = {}, set()
    for g in todo:
        c = select_candidates(reps, pf, groups[g], cfg)
        tabs[g] = c
        union_ids |= set(c["protein_id"].dropna())
    if todo:
        allfa = scratch / "cand_all.faa"
        n = extract_from_db(stage / "genus_reps_db" / "genus_reps", union_ids, allfa, scratch / "ex")
        print(f"extracted {n} of {len(union_ids)} candidate sequences", flush=True)
        for g in todo:
            c = tabs[g]
            ids = set(c["protein_id"].dropna())
            nn = filter_fasta(allfa, ids, stage / f"cand_{g}.faa.gz")
            c.to_csv(stage / f"cand_{g}.tsv.gz", sep="\t", index=False)
            print(f"candidates {g} ({groups[g]}): {len(ids)} proteins in {int(c.loc[c.present, 'species_prefix'].nunique())} genera, "
                  f"{int((~c.present).sum())} genera without the domain; {nn} sequences", flush=True)
        allfa.unlink()


def step_search(stage: Path, root: Path, scratch: Path, min_pident=MIN_PIDENT, min_cov=MIN_COV, cov_mode=COV_MODE):
    cfg = load_cfg(root)
    sdf = pd.read_csv(stage / "seeds.tsv", sep="\t", keep_default_na=False)
    counts = []
    for grp in sorted(set(GROUP.values())):
        tg = [t for t in GROUP if GROUP[t] == grp]
        hp, ap = stage / f"hits_{grp}.tsv", stage / f"assigned_{grp}.tsv"
        cand = pd.read_csv(stage / f"cand_{grp}.tsv.gz", sep="\t", keep_default_na=True)
        cand["present"] = cand["present"].astype(bool)
        rules = {t: sd.rules_for_target(cfg, t)[1] for t in tg}
        params = {"min_pident": min_pident, "min_cov": min_cov, "cov_mode": cov_mode}
        if ap.exists() and hp.exists():
            hits = pd.read_csv(hp, sep="\t", comment="#")
            asg = pd.read_csv(ap, sep="\t", comment="#")
        else:
            hits = search(stage / f"seeds_{grp}.faa", stage / f"cand_{grp}.faa.gz", hp, scratch / grp, min_pident, min_cov,
                          cov_mode=cov_mode)
            st = dict(zip(sdf["seed_id"], sdf["target"]))
            asg = assign_best_seed(hits, cand, st, rules)
            write_commented(asg, ap, dict(params, rules="; ".join(f"{t}: {'&'.join(r) or 'none'}" for t, r in rules.items())))
        for t in tg:
            slug = SLUG[t]
            ph = pass_for_target(asg, t)
            write_commented(ph, stage / f"pass_{slug}.tsv", dict(params, architecture_rule="&".join(rules[t]) or "none"))
            sub = sdf[sdf["target"] == t]
            mine = asg[asg["seed_target"] == t]
            counts.append({"target": t, "n_seeds": len(sub), "n_seeds_ok": int((sub.status == "ok").sum()),
                           "n_candidate_proteins": int(cand["protein_id"].notna().sum()),
                           "n_candidate_genera": int(cand.loc[cand.present, "species_prefix"].nunique()),
                           "n_reps": int(cand["species_prefix"].nunique()),
                           "n_hit_proteins": int(hits["target"].nunique()) if len(hits) else 0,
                           "n_best_seed_of_target": len(mine),
                           "n_arch_pass_proteins": len(ph), "n_arch_pass_genera": int(ph["species_prefix"].nunique()) if len(ph) else 0,
                           "n_disagree": count_disagreements(asg, t),
                           "n_disagree_other_type": int(((mine["rule_ok"] == False) & (mine["arch_type"] != "")).sum()) if len(mine) else 0,  # noqa: E712
                           "min_pident": min_pident, "min_cov": min_cov, "cov_mode": cov_mode,
                           "architecture_rule": "&".join(rules[t]) or "none"})
            print(counts[-1], flush=True)
    pd.DataFrame(counts).to_csv(stage / "counts.tsv", sep="\t", index=False)


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("step", choices=["seeds", "candidates", "search"])
    ap.add_argument("--stage", default=str(RUNS / "analysis/tf_orthologs"))
    ap.add_argument("--targets-stage", default=str(RUNS / "analysis/tf_targets"), help="holds genus_reps.tsv, genus_reps_db, FBF4BF29.faa")
    ap.add_argument("--tables", default=str(RUNS / "tables"), help="BFD parquet tables (read only)")
    ap.add_argument("--root", default=str(HERE.parents[1]))
    ap.add_argument("--scratch", required=True)
    ap.add_argument("--threads", type=int, default=None)
    a = ap.parse_args(argv)
    stage, root, scratch = Path(a.stage), Path(a.root), Path(a.scratch)
    scratch.mkdir(parents=True, exist_ok=True)
    stage.mkdir(parents=True, exist_ok=True)
    tst = Path(a.targets_stage)
    for f in ("genus_reps.tsv", "FBF4BF29.faa", "genus_reps_db"):      # reuse the Task 5b files: link, do not rebuild
        link_reuse(stage / f, tst / f)
    if invalidate_if_changed(stage, root):
        print("seed config, thresholds or code changed (or first run): derived outputs removed, steps rerun", flush=True)
    if a.step == "seeds":
        step_seeds(stage, root, Path(a.tables), scratch)
    elif a.step == "candidates":
        step_candidates(stage, root, Path(a.tables), scratch)
    else:
        step_search(stage, root, scratch)


if __name__ == "__main__":
    main(sys.argv[1:])
