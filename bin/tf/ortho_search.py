"""Ortholog search for WC-2, Zap1-like and bZip5 among the genus-representative proteomes.

Steps (each skips when its output exists; run from the repo root):
    python bin/tf/ortho_search.py seeds       seed table, seed FASTA per target, best BFD match of every seed
    python bin/tf/ortho_search.py candidates  proteins of the genus representatives with the target Pfam, with architecture
    python bin/tf/ortho_search.py search      MMseqs2 seeds vs candidates, then the architecture rule

Order of the rules: candidates carry the target Pfam only. The search (identity and coverage thresholds) runs on all of
them. The architecture rule (for example GATA plus PAS) is applied to the hits AFTER the search. A genus
representative without the target domain is kept in the candidate table as an absent row, not dropped.
Similarity and tree position are evidence for orthology, not proof.
"""
import argparse
import gzip
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
SLUG = {"WC-2": "wc2", "Zap1-like": "zap1", "bZip5": "bzip5"}
ANCHORS = {"FBF4BF29": "Rhizopus microsporus ATCC 11559", "F750FE6D": "Mucor lusitanicus", "F43EADC6": "Phycomyces blakesleeanus"}
BFD_GENOME = "FBF4BF29"
HIT_COLS = jm.HIT_COLS


# ---------------------------------------------------------------- rules
def architecture_ok(arch: str, required) -> bool:
    """True when the architecture string (accessions joined by '>', no versions) meets every entry of `required`.
    An entry 'A|B' needs at least one domain of A or B. An entry 'A>=N' needs at least N domains of A.
    An empty list is always true."""
    doms = [d for d in str(arch).split(">") if d] if isinstance(arch, str) else []
    for entry in required:
        if ">=" in entry:
            acc, n = entry.split(">=")
            if sum(1 for d in doms if d in set(acc.split("|"))) < int(n):
                return False
        elif not (set(entry.split("|")) & set(doms)):
            return False
    return True


def architecture_names(arch, names):
    """Accession architecture -> readable text with Pfam names (names: accession -> name)."""
    return ">".join(names.get(a, a) for a in str(arch).split(">") if a)


# ---------------------------------------------------------------- candidates
def select_candidates(genus_reps: pd.DataFrame, pfam_hits: pd.DataFrame, target: str, seeds: pd.DataFrame) -> pd.DataFrame:
    """All proteins of the genus representatives that carry the target Pfam (seeds['target_pfam'] of `target`).
    pfam_hits: every Pfam domain row of those proteins (protein_id, pfam_acc, ali_from). Columns of the result:
    species_prefix, protein_id, architecture, n_target_domains, present. A representative with no such protein
    gets one row with protein_id NaN, empty architecture and present False (reported as absent)."""
    tpf = seeds.loc[seeds["target"] == target, "target_pfam"].iloc[0]
    ph = pfam_hits.copy()
    ph["acc"] = ph["pfam_acc"].map(sd.strip_version)
    ph = ph.sort_values(["protein_id", "ali_from"], kind="stable")
    arch = ph.groupby("protein_id")["acc"].agg(">".join)
    ntar = ph[ph["acc"] == tpf].groupby("protein_id").size()
    reps = set(genus_reps["species_prefix"])
    rows = [{"species_prefix": str(p).split("_")[0], "protein_id": p, "architecture": arch[p],
             "n_target_domains": int(n), "present": True} for p, n in ntar.items() if str(p).split("_")[0] in reps]
    df = pd.DataFrame(rows, columns=["species_prefix", "protein_id", "architecture", "n_target_domains", "present"])
    have = set(df["species_prefix"])
    absent = pd.DataFrame([{"species_prefix": g, "protein_id": np.nan, "architecture": "", "n_target_domains": 0, "present": False}
                           for g in genus_reps["species_prefix"] if g not in have],
                          columns=df.columns)
    out = pd.concat([df, absent], ignore_index=True)
    return out.sort_values(["species_prefix", "protein_id"], na_position="last", kind="stable").reset_index(drop=True)


# ---------------------------------------------------------------- search
def write_commented(df: pd.DataFrame, path, params: dict):
    """Write a TSV whose first line is '# key=value ...' with the thresholds. Read it with comment='#'."""
    with open(path, "w") as f:
        f.write("# " + " ".join(f"{k}={v}" for k, v in params.items()) + "\n")
        df.to_csv(f, sep="\t", index=False)


def search(seed_fasta: Path, cand_fasta: Path, out: Path, scratch: Path, min_pident=MIN_PIDENT, min_cov=MIN_COV,
           threads=None) -> pd.DataFrame:
    """MMseqs2 easy-search of the seeds against the candidates. Thresholds: identity min_pident (percent), coverage
    min_cov of both sequences (cov mode 0). The thresholds are written as a '# ' line at the top of `out` and kept in
    df.attrs. Without mmseqs on PATH the step is skipped: an empty table, attrs['skipped'] says why."""
    params = {"min_pident": float(min_pident), "min_cov": float(min_cov), "cov_mode": 0, "max_seqs": 100000}
    if shutil.which("mmseqs") is None:
        msg = "mmseqs is not on PATH: search skipped (module load MMseqs2/17-b804f)"
        print(msg, flush=True)
        df = pd.DataFrame(columns=HIT_COLS)
        df.attrs.update(params, skipped=msg)
        return df
    Path(scratch).mkdir(parents=True, exist_ok=True)
    raw = Path(str(out) + ".raw")
    df = jm.easy_search(Path(seed_fasta), Path(cand_fasta), raw, Path(scratch), min_pident, min_cov, threads=threads,
                        max_seqs=params["max_seqs"])
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
    m = best.merge(cand[cand["present"]][["species_prefix", "protein_id", "architecture"]],
                   left_on="target", right_on="protein_id", how="inner")
    m = m[m["architecture"].map(lambda a: architecture_ok(a, required))]
    m = m.rename(columns={"query": "seed"})
    return m[cols].sort_values(["species_prefix", "protein_id"]).reset_index(drop=True)


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


def step_seeds(stage: Path, root: Path, tables_dir: Path, scratch: Path):
    out = stage / "seeds.tsv"
    if out.exists():
        print(f"seeds: {out} exists, skipped", flush=True)
        return pd.read_csv(out, sep="\t", keep_default_na=False)
    cfg = load_cfg(root)
    paths = None
    rows = []
    # paper seeds (JGI full-length sequences)
    pap = sd.paper_sequences(cfg, JGI_FASTA)
    for (_, c), (_, p) in zip(cfg[cfg["source"] == "paper"].iterrows(), pap.iterrows()):
        rows.append({"target": c["target"], "label": c["seed_name"], "source": "paper", "key": c["key"], "seed_id": p["accession"],
                     "accession": p["accession"], "sequence": p["sequence"], "status": p["status"], "genome": BFD_GENOME,
                     "phylum": "Mucoromycota"})
    # UniProt seeds, cached
    uni = sd.fetch_seeds(cfg, stage / "seed_cache")
    for (_, c), (_, u) in zip(cfg[cfg["source"] == "uniprot"].iterrows(), uni.iterrows()):
        rows.append({"target": c["target"], "label": c["seed_name"], "source": "uniprot", "key": c["key"],
                     "seed_id": u["accession"] or f"uniprot_{c['key']}", "accession": u["accession"], "sequence": u["sequence"],
                     "status": u["status"], "genome": "", "phylum": sd.TAXON_PHYLUM.get(c["taxon_id"], "")})
    # BFD seeds selected by domain
    bfd_rows = cfg[cfg["source"] == "bfd"]
    if len(bfd_rows):
        paths = jm_proteome_paths()
        gen = [k.split("|")[0] for k in bfd_rows["key"]]
        pf = read_pfam_for_genomes(tables_dir, gen)
        for _, c in bfd_rows.iterrows():
            g, sel = c["key"].split("|")
            ids = sd.bfd_selector_proteins(g, sel, pf)
            seqs = {}
            if g in paths and Path(paths[g]).exists():
                seqs = sd.read_fasta_dict(paths[g])
            if not ids:
                rows.append({"target": c["target"], "label": c["seed_name"], "source": "bfd", "key": c["key"], "seed_id": "",
                             "accession": "", "sequence": "", "status": f"missing: no protein of {g} matches {sel}",
                             "genome": g, "phylum": "Mucoromycota"})
            for pid in ids:
                rows.append({"target": c["target"], "label": f"{c['seed_name']} {pid}", "source": "bfd", "key": c["key"],
                             "seed_id": pid, "accession": pid, "sequence": seqs.get(pid, ""),
                             "status": "ok" if seqs.get(pid) else "missing: protein sequence not found", "genome": g,
                             "phylum": "Mucoromycota"})
    sdf = pd.DataFrame(rows)
    # best BFD (FBF4BF29) match of every seed with the thresholds of the ortholog search, and its architecture
    bfd_faa = stage / "FBF4BF29.faa"
    ok = sdf[sdf["status"] == "ok"]
    q = scratch / "all_seeds.faa"
    q.write_text("".join(f">{r.seed_id}\n{r.sequence}\n" for r in ok.itertuples()))
    hits = jm.easy_search(q, bfd_faa, scratch / "seeds_vs_bfd.tsv", scratch / "mm", MIN_PIDENT, MIN_COV, max_seqs=20)
    best = hits.sort_values("bits", ascending=False).drop_duplicates("query").set_index("query") if len(hits) else pd.DataFrame()
    pfb = read_pfam_for_genomes(tables_dir, [BFD_GENOME])
    arch = arch_map(pfb)
    names = pd.read_parquet(RUNS / "analysis/summary_tables/pfam_names.parquet").set_index("pfam_acc")["pfam_id"].to_dict()
    mt, mp, mq, ma, man = [], [], [], [], []
    for r in sdf.itertuples():
        if r.seed_id in best.index:
            b = best.loc[r.seed_id]
            a = arch.get(b["target"], "")
            mt.append(b["target"]); mp.append(float(b["pident"])); mq.append(float(min(b["qcov"], b["tcov"]))); ma.append(a)
            man.append(architecture_names(a, names))
        else:
            mt.append(""); mp.append(np.nan); mq.append(np.nan); ma.append(""); man.append("")
    sdf["bfd_match"], sdf["bfd_pident"], sdf["bfd_min_cov"] = mt, mp, mq
    sdf["bfd_architecture"], sdf["bfd_architecture_names"] = ma, man
    # paper-named sanity check: resolve against BFD at 95/0.9 with the paper names
    pap_seqs = {r.seed_id: r.sequence for r in sdf[(sdf.source == "paper") & (sdf.status == "ok")].itertuples()}
    res = sd.resolve_paper_seeds(pap_seqs, bfd_faa, ROOT_NAMES(root), arch, scratch / "resolve",
                               min_pident=90.0, min_cov=0.8)   # as the Task 5b JGI-to-BFD map: gene models differ in length
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
            rule = [tpf] + req if tpf else req
            if not r["bfd_match"]:
                w = f"no BFD {BFD_GENOME} match at identity {MIN_PIDENT} and coverage {MIN_COV}"
            elif not architecture_ok(r["bfd_architecture"], rule):
                w = (f"architecture check failed: best BFD match {r['bfd_match']} has '{r['bfd_architecture']}', "
                     f"rule for {r['target']}: {rule}")
        sdf.loc[i, "warning"] = w
    sdf.to_csv(out, sep="\t", index=False)
    for t, g in sdf[sdf["status"] == "ok"].groupby("target"):
        (stage / f"seeds_{SLUG[t]}.faa").write_text("".join(f">{r.seed_id}\n{r.sequence}\n" for r in g.itertuples()))
    for t in SLUG:
        sub = sdf[sdf["target"] == t]
        print(f"seeds {t}: {int((sub.status == 'ok').sum())} ok, {int((sub.status != 'ok').sum())} missing", flush=True)
    return sdf


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
    return con.execute(f"select protein_id, pfam_id, pfam_acc, ali_from from {pq} where split_part(protein_id, '_', 1) in ({gl})").df()


def step_candidates(stage: Path, root: Path, tables_dir: Path, scratch: Path):
    cfg = load_cfg(root)
    reps = pd.read_csv(stage / "genus_reps.tsv", sep="\t")
    accs = sorted({sd.strip_version(a) for t in SLUG for a in [sd.rules_for_target(cfg, t)[0]]})
    pf_path = stage / "pfam_candidates.parquet"
    if pf_path.exists():
        pf = pd.read_parquet(pf_path)
    else:
        pf = read_pfam_rows(tables_dir, reps["species_prefix"], accs)
        pf.to_parquet(pf_path, index=False)
    print(f"pfam rows of candidate proteins: {len(pf)}", flush=True)
    todo = [t for t in SLUG if not (stage / f"cand_{SLUG[t]}.faa.gz").exists()]
    union_ids, tabs = set(), {}
    for t in todo:
        c = select_candidates(reps, pf, t, cfg)
        tabs[t] = c
        union_ids |= set(c["protein_id"].dropna())
    if todo:
        allfa = scratch / "cand_all.faa"
        n = extract_from_db(stage / "genus_reps_db" / "genus_reps", union_ids, allfa, scratch / "ex")
        print(f"extracted {n} of {len(union_ids)} candidate sequences", flush=True)
        for t in todo:
            c = tabs[t]
            ids = set(c["protein_id"].dropna())
            nn = filter_fasta(allfa, ids, stage / f"cand_{SLUG[t]}.faa.gz")
            c.to_csv(stage / f"cand_{SLUG[t]}.tsv.gz", sep="\t", index=False)
            print(f"candidates {t}: {len(ids)} proteins in {int(c.loc[c.present, 'species_prefix'].nunique())} genera, "
                  f"{int((~c.present).sum())} genera without the domain; {nn} sequences", flush=True)
        allfa.unlink()


def step_search(stage: Path, root: Path, scratch: Path, min_pident=MIN_PIDENT, min_cov=MIN_COV):
    cfg = load_cfg(root)
    sdf = pd.read_csv(stage / "seeds.tsv", sep="\t", keep_default_na=False)
    counts = []
    for t, slug in SLUG.items():
        hp = stage / f"hits_{slug}.tsv"
        pp = stage / f"pass_{slug}.tsv"
        cand = pd.read_csv(stage / f"cand_{slug}.tsv.gz", sep="\t", keep_default_na=True)
        cand["present"] = cand["present"].astype(bool)
        _, req = sd.rules_for_target(cfg, t)
        if pp.exists():
            ph = pd.read_csv(pp, sep="\t", comment="#")
            hits = pd.read_csv(hp, sep="\t", comment="#") if hp.exists() else pd.DataFrame()
        else:
            hits = search(stage / f"seeds_{slug}.faa", stage / f"cand_{slug}.faa.gz", hp, scratch / slug, min_pident, min_cov)
            ph = pass_hits(hits, cand, req)
            write_commented(ph, pp, {"min_pident": min_pident, "min_cov": min_cov, "architecture_rule": ";".join(req) or "none"})
        sub = sdf[sdf["target"] == t]
        counts.append({"target": t, "n_seeds": len(sub), "n_seeds_ok": int((sub.status == "ok").sum()),
                       "n_candidate_proteins": int(cand["protein_id"].notna().sum()),
                       "n_candidate_genera": int(cand.loc[cand.present, "species_prefix"].nunique()),
                       "n_reps": int(cand["species_prefix"].nunique()),
                       "n_hit_proteins": int(hits["target"].nunique()) if len(hits) else 0,
                       "n_hit_genera": int(cand[cand.protein_id.isin(set(hits["target"]))]["species_prefix"].nunique()) if len(hits) else 0,
                       "n_arch_pass_proteins": len(ph), "n_arch_pass_genera": int(ph["species_prefix"].nunique()) if len(ph) else 0,
                       "min_pident": min_pident, "min_cov": min_cov, "architecture_rule": ";".join(req) or "none"})
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
    for f in ("genus_reps.tsv", "FBF4BF29.faa"):      # reuse the Task 5b files: link, do not rebuild
        if not (stage / f).exists():
            (stage / f).symlink_to(tst / f)
    if not (stage / "genus_reps_db").exists():
        (stage / "genus_reps_db").symlink_to(tst / "genus_reps_db")
    if a.step == "seeds":
        step_seeds(stage, root, Path(a.tables), scratch)
    elif a.step == "candidates":
        step_candidates(stage, root, Path(a.tables), scratch)
    else:
        step_search(stage, root, scratch)


if __name__ == "__main__":
    main(sys.argv[1:])
