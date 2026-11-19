"""Annotation-gap check: does a GATA+PAS gene exist outside the gene models of a genome?

Translated search (MMseqs2, protein seeds against the genome assembly, --search-type 2) of the WC-2 and WC-1 seeds. Each
hit is compared with the gene models (GFF3) of the genome: it overlaps a model (the model id and its Pfam architecture
are reported) or it lies in no model. This is evidence only. It does not show that a hit is a functional gene.

    python bin/tf/ortho_gap_check.py --stage DIR --genomes F07749BD F8122E67 F75D53DE --scratch DIR
"""
import argparse
import gzip
import subprocess
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import jgi_map as jm  # noqa: E402
import ortho_search as ors  # noqa: E402
import seeds as sd  # noqa: E402

RUNS = ors.RUNS
FMT = "query,target,pident,alnlen,qcov,qstart,qend,tstart,tend,evalue,bits"
COLS = FMT.split(",")


def read_genes(gff3_gz):
    """Gene models of a funannotate GFF3: DataFrame scaffold, start, end, gene_id, protein_id (mRNA id)."""
    rows, mr = [], {}
    with gzip.open(gff3_gz, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            at = dict(x.split("=", 1) for x in c[8].split(";") if "=" in x)
            if c[2] == "gene":
                rows.append({"scaffold": c[0], "start": int(c[3]), "end": int(c[4]), "gene_id": at.get("ID", "")})
            elif c[2] == "mRNA":
                mr.setdefault(at.get("Parent", ""), at.get("ID", ""))
    df = pd.DataFrame(rows, columns=["scaffold", "start", "end", "gene_id"])
    df["protein_id"] = df["gene_id"].map(mr)
    return df


def overlap_genes(hits: pd.DataFrame, genes: pd.DataFrame) -> pd.DataFrame:
    """Add gene_id and protein_id of the gene model that overlaps each hit (tstart/tend on the target scaffold; the
    larger overlap wins). Hits in no model get empty ids."""
    out = hits.copy()
    gid, pid = [], []
    by = {k: g for k, g in genes.groupby("scaffold")}
    for r in out.itertuples():
        lo, hi = sorted((int(r.tstart), int(r.tend)))
        g = by.get(r.target)
        best, bo = None, 0
        if g is not None:
            for x in g.itertuples():
                o = min(hi, x.end) - max(lo, x.start) + 1
                if o > bo:
                    best, bo = x, o
        gid.append(best.gene_id if best is not None else "")
        pid.append(best.protein_id if best is not None else "")
    out["gene_id"], out["protein_id"] = gid, pid
    return out


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", required=True)
    ap.add_argument("--genomes", nargs="+", required=True)
    ap.add_argument("--scratch", required=True)
    ap.add_argument("--tables", default=str(RUNS / "tables"))
    ap.add_argument("--targets", nargs="+", default=["WC-2", "WC-1"])
    a = ap.parse_args(argv)
    stage, scratch = Path(a.stage), Path(a.scratch)
    scratch.mkdir(parents=True, exist_ok=True)
    seeds = pd.read_csv(stage / "seeds.tsv", sep="\t", keep_default_na=False)
    keep = seeds[seeds["target"].isin(a.targets) & (seeds["status"] == "ok") & seeds["source"].isin(["paper", "uniprot"])]
    q = scratch / "gap_seeds.faa"
    q.write_text("".join(f">{r.seed_id}|{r.target}\n{r.sequence}\n" for r in keep.itertuples()))
    names = {g: n for g, n in (l.split("\t")[:2] for l in open(RUNS / "seqdb/genomes.tsv"))}
    exe = jm.need_mmseqs()
    allrows = []
    for g in a.genomes:
        name = names[g]
        pr = RUNS / "genome_annotation" / name / "predict_results"
        fa, gff = pr / f"{name}.scaffolds.fa.gz", pr / f"{name}.gff3.gz"
        raw = scratch / f"{g}.tsv"
        subprocess.run([exe, "easy-search", str(q), str(fa), str(raw), str(scratch / f"mm_{g}"), "--search-type", "2",
                        "-e", "1e-5", "--format-output", FMT, "--threads", "4", "-v", "1"], check=True)
        h = pd.read_csv(raw, sep="\t", header=None, names=COLS) if raw.stat().st_size else pd.DataFrame(columns=COLS)
        h = overlap_genes(h, read_genes(gff))
        pf = ors.read_pfam_for_genomes(a.tables, [g])
        arch, npas = ors.arch_map(pf), ors.npas_map(pf)
        h["architecture"] = h["protein_id"].map(lambda p: arch.get(p, ""))
        h["n_pas_regions"] = h["protein_id"].map(lambda p: npas.get(p, 0))
        h.insert(0, "genome", g)
        allrows.append(h)
    res = pd.concat(allrows, ignore_index=True)
    res.to_csv(stage / "gap_check_hits.tsv", sep="\t", index=False)
    # summary: best hit per genome and seed, and whether a hit lies in a model with GATA and PAS
    top = res.sort_values("bits", ascending=False).drop_duplicates(["genome", "query"])
    print(top.to_string(index=False))
    gata = res[res["architecture"].str.contains("PF00320", na=False)]
    print("\nhits in a gene model with GATA (PF00320):", len(gata), "; models:", sorted(set(gata["protein_id"])))
    print("hits in no gene model:", int((res["gene_id"] == "").sum()), "of", len(res))


if __name__ == "__main__":
    main(sys.argv[1:])
