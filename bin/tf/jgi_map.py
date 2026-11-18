"""JGI (Rhimi59_2) protein readers and the JGI-to-BFD protein ID map for R. microsporus ATCC 11559.

The DAP-seq targets of Lax et al. 2025 are JGI protein IDs. BFD has its own annotation (genome FBF4BF29,
assembly GCA_002083735.1), so the IDs are mapped by reciprocal best hits (RBH) with MMseqs2.
"""
import gzip
import os
import re
import shutil
import subprocess
from pathlib import Path

import openpyxl
import pandas as pd

MIN_PIDENT = 90.0
MIN_COV = 0.8
FORMAT = "query,target,pident,qcov,tcov,evalue,bits"
HIT_COLS = FORMAT.split(",")


def _open_text(p):
    p = str(p)
    if p.endswith(".gz"):
        return gzip.open(p, "rt")
    if p.endswith(".zst"):
        return subprocess.Popen(["zstd", "-dc", p], stdout=subprocess.PIPE, text=True).stdout
    return open(p)


# ---------------------------------------------------------------- readers
def read_jgi_proteins(fasta_gz: Path) -> dict:
    """proteinId -> (name, sequence). Header: jgi|<portal>|<proteinId>|<name>."""
    out, pid, name, seq = {}, None, "", []
    with _open_text(fasta_gz) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if pid is not None:
                    out[pid] = (name, "".join(seq))
                parts = line[1:].split("|", 3)
                pid = parts[2] if len(parts) >= 3 else parts[0]
                name = parts[3] if len(parts) > 3 else ""
                seq = []
            elif line:
                seq.append(line.strip())
    if pid is not None:
        out[pid] = (name, "".join(seq))
    return out


def _attrs(text):
    d = {}
    for kv in text.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = v
    return d


def read_jgi_catalog(gff3_gz: Path) -> pd.DataFrame:
    """One row per gene of the filtered best models. track comes from the mRNA line (FilteredModelsN)."""
    genes, track = {}, {}
    with _open_text(gff3_gz) as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            a = _attrs(c[8])
            if c[2] == "gene":
                pid = a.get("proteinId")
                if pid is None:
                    continue
                genes[pid] = (pid, c[0], int(c[3]), int(c[4]), c[6], a.get("product_name", ""))
            elif c[2] == "mRNA" and "proteinId" in a:
                track[a["proteinId"]] = a.get("track", "")
    rows = [g + (track.get(pid, ""),) for pid, g in genes.items()]
    return pd.DataFrame(rows, columns=["protein_id", "scaffold", "start", "end", "strand", "product_name", "track"])


def _norm(name):
    return re.sub(r"[\s_]+", "", str(name)).lower()


def read_tf_names(tsv: Path) -> pd.DataFrame:
    """config/paper_tf_names.tsv: 88 tested TFs. IDs stay strings."""
    return pd.read_csv(tsv, sep="\t", comment="#", dtype={"jgi_id": str})


def read_targets(xlsx: Path, names: pd.DataFrame = None) -> pd.DataFrame:
    """One row per (TF sheet, target protein). Sheet names are kept as written (stripped). Duplicates within a
    sheet are removed. With a names table (read_tf_names) each row also gets tf_jgi_id and tf_family; a sheet
    that matches no name (spaces and underscores removed, case ignored) raises ValueError."""
    wb = openpyxl.load_workbook(xlsx, read_only=True, data_only=True)
    rows = []
    try:
        for sheet in wb.sheetnames:
            ws = wb[sheet]
            it = ws.iter_rows(values_only=True)
            header = next(it, None)
            if header is None:
                continue
            col = next((i for i, h in enumerate(header) if h and str(h).strip().startswith("Target")), 0)
            seen = set()
            for r in it:
                v = r[col] if col < len(r) else None
                if v is None or str(v).strip() == "":
                    continue
                pid = str(int(v)) if isinstance(v, float) and v == int(v) else str(v).strip()
                if pid not in seen:
                    seen.add(pid)
                    rows.append((sheet.strip(), pid))
    finally:
        wb.close()
    df = pd.DataFrame(rows, columns=["tf_name", "jgi_protein_id"])
    if names is not None:
        norm = [_norm(n) for n in names["name"]]

        def find(sheet):
            k = _norm(sheet)
            hit = [i for i, n in enumerate(norm) if n == k] or [i for i, n in enumerate(norm) if n.startswith(k)]
            return hit[0] if len(hit) == 1 else None    # exact match first; else a unique name that starts with it (CP2 -> CP2_TF1)
        pos = {s: find(s) for s in df["tf_name"].unique()}
        bad = sorted(s for s, i in pos.items() if i is None)
        if bad:
            raise ValueError(f"sheet names without a match in the TF names table: {bad}")
        idx = df["tf_name"].map(pos)
        df["tf_jgi_id"] = names["jgi_id"].astype(str).to_numpy()[idx]
        df["tf_family"] = names["tf_family"].to_numpy()[idx]
    return df


def read_tested_tfs(xlsx: Path) -> pd.DataFrame:
    """Tested TF fragments (MOESM6): Rhimi59_2_<id>[a|b] -> base protein ID and fragment suffix."""
    wb = openpyxl.load_workbook(xlsx, read_only=True, data_only=True)
    rows = []
    try:
        for ws in wb.worksheets:
            for r in ws.iter_rows(values_only=True):
                v = r[0] if r else None
                m = re.fullmatch(r"Rhimi59_2_(\d+)([a-z]?)", str(v).strip()) if v else None
                if m:
                    rows.append((str(v).strip(), m.group(1), m.group(2)))
    finally:
        wb.close()
    return pd.DataFrame(rows, columns=["tested_id", "jgi_protein_id", "fragment"])


# ---------------------------------------------------------------- MMseqs2 wrappers
def need_mmseqs():
    exe = shutil.which("mmseqs")
    if exe is None:
        raise RuntimeError("mmseqs is not on PATH. Run: module load MMseqs2/17-b804f")
    return exe


def threads_default():
    return int(os.environ.get("SLURM_CPUS_PER_TASK", "4"))


def easy_search(query: Path, target: Path, out_tsv: Path, tmp: Path, min_pident, min_cov, threads=None,
                split_memory_limit=None, max_seqs=50):
    """MMseqs2 easy-search; writes a table with the columns in FORMAT. pident is 0-100, qcov/tcov are 0-1."""
    exe = need_mmseqs()
    tmp.mkdir(parents=True, exist_ok=True)
    cmd = [exe, "easy-search", str(query), str(target), str(out_tsv), str(tmp / "mm_tmp"),
           "--min-seq-id", str(min_pident / 100.0), "-c", str(min_cov), "--cov-mode", "0",
           "--format-output", FORMAT, "--threads", str(threads or threads_default()), "--max-seqs", str(max_seqs),
           "-v", "1"]
    if split_memory_limit:
        cmd += ["--split-memory-limit", str(split_memory_limit)]
    subprocess.run(cmd, check=True)
    shutil.rmtree(tmp / "mm_tmp", ignore_errors=True)
    return read_hits(out_tsv)


def read_hits(path):
    if Path(path).stat().st_size == 0:
        return pd.DataFrame(columns=HIT_COLS)
    df = pd.read_csv(path, sep="\t", header=None, names=HIT_COLS, dtype={"query": str, "target": str})
    return df


def fasta_ids(path):
    with _open_text(path) as f:
        return [line[1:].split()[0] for line in f if line.startswith(">")]


def rbh_table(ids, fwd: pd.DataFrame, rev: pd.DataFrame) -> pd.DataFrame:
    """Reciprocal best hits from a forward (JGI -> BFD) and reverse (BFD -> JGI) hit table. The best hit is the
    highest bit score. Hits tied at the top bit score all count as best (identical sequences). A JGI protein is
    mapped when its best BFD hit has the JGI protein among the best hits of the reverse search; among tied
    forward hits the first reciprocal one is taken. Others get an empty bfd_protein_id and rbh False."""
    rev_best = {}
    for q, g in rev.groupby("query"):
        top = g["bits"].max()
        rev_best[q] = set(g.loc[g["bits"] >= top, "target"])
    fwd_best = {}
    for q, g in fwd.groupby("query"):
        top = g["bits"].max()
        fwd_best[q] = g[g["bits"] >= top].sort_values(["bits", "target"], ascending=[False, True])
    rows = []
    for q in ids:
        hit = None
        for _, r in fwd_best.get(q, pd.DataFrame()).iterrows():
            if q in rev_best.get(r["target"], ()):
                hit = r
                break
        if hit is None:
            rows.append((q, "", float("nan"), float("nan"), float("nan"), False))
        else:
            rows.append((q, hit["target"], float(hit["pident"]), float(hit["qcov"]), float(hit["tcov"]), True))
    return pd.DataFrame(rows, columns=["jgi_protein_id", "bfd_protein_id", "pident", "qcov", "tcov", "rbh"])


def map_jgi_to_bfd(jgi_fasta: Path, bfd_fasta: Path, out_tsv: Path, scratch: Path, min_pident=MIN_PIDENT,
                   min_cov=MIN_COV, threads=None) -> pd.DataFrame:
    """MMseqs2 easy-search both ways and keep reciprocal best hits. Thresholds go to the first line of out_tsv.
    All JGI proteins of the FASTA are in the result; those without an RBH have an empty bfd_protein_id."""
    scratch = Path(scratch)
    f = easy_search(jgi_fasta, bfd_fasta, scratch / "fwd.tsv", scratch / "fwd", min_pident, min_cov, threads)
    r = easy_search(bfd_fasta, jgi_fasta, scratch / "rev.tsv", scratch / "rev", min_pident, min_cov, threads)
    df = rbh_table(fasta_ids(jgi_fasta), f, r)
    out_tsv = Path(out_tsv)
    with open(out_tsv, "w") as fh:
        fh.write(f"# min_pident={min_pident} min_cov={min_cov} cov_mode=0 (both sequences) rbh=best bit score, ties allowed\n")
        df.to_csv(fh, sep="\t", index=False)
    return df


def read_map(path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", comment="#", dtype={"jgi_protein_id": str, "bfd_protein_id": str},
                     keep_default_na=False, na_values=[""])
    df["bfd_protein_id"] = df["bfd_protein_id"].fillna("")
    return df


def mapping_rate(map_df: pd.DataFrame, ids, best_ids=None) -> dict:
    """Share of ids with a BFD partner. best_ids: the filtered best models of the GFF3; the same numbers are
    given for the ids inside it (n_ids_best, n_mapped_best, rate_best)."""
    mapped = set(map_df.loc[map_df["bfd_protein_id"].fillna("") != "", "jgi_protein_id"].astype(str))
    ids = list(dict.fromkeys(str(i) for i in ids))
    n = len(ids)
    k = sum(i in mapped for i in ids)
    out = {"n_ids": n, "n_mapped": k, "rate": k / n if n else float("nan")}
    if best_ids is not None:
        bset = set(map(str, best_ids))
        b = [i for i in ids if i in bset]
        kb = sum(i in mapped for i in b)
        out.update({"n_ids_best": len(b), "n_mapped_best": kb, "rate_best": kb / len(b) if b else float("nan")})
    return out
