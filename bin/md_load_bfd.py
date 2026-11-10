#!/usr/bin/env python3
"""Load BFD parquet tables into a MotherDuck database.

Each parquet file (or directory of parquet parts) in --tables becomes one table.
Table names: file stem with '.' -> '_' (cazy.overview -> cazy_overview).

Restartable: a table that already exists in MotherDuck with the same row count
as its source is skipped, so a failed or timed-out job can simply be resubmitted.

Table comments come from --tables/_manifest/<name>.json ("description" key)
when present, else from a short default.

Token: read from $motherduck_token, else from --token-file
(default ~/.config/motherduck/token, must be chmod 600). Never pass it on the
command line.

Usage:
    bin/md_load_bfd.py --tables /bigdata/.../Fungi_BFD_runs/tables --db bfd --dry-run
    bin/md_load_bfd.py --tables /bigdata/.../Fungi_BFD_runs/tables --db bfd
    bin/md_load_bfd.py ... --only pfam species        # subset
    bin/md_load_bfd.py ... --include-swissprot         # also load swissprot_hits (~37 GB)
"""
import argparse
import fnmatch
import json
import os
import stat
import sys
import time
from pathlib import Path

import duckdb

# Excluded unless --include-swissprot / explicitly listed with --only.
DEFAULT_EXCLUDE = ["swissprot_hits", "*.bak*", "*.csv*", "*.md", "_*"]

DEFAULT_COMMENTS = {
    "species": "One row per genome. Taxonomy PHYLUM..SPECIES. Keys: ASMID; LOCUSTAG = protein_id prefix = species_prefix.",
    "samples": "Per-genome sample metadata (BioProject, strain, translation table). Key ASMID.",
    "asm_stats": "Assembly statistics per genome (ASMID): N50/L50, GC, repeat masking, telomeres, annotation status.",
    "busco_genome": "BUSCO completeness per genome (ASMID).",
    "pfam": "Pfam hmmsearch domain hits, one row per domain. Genome = split_part(protein_id,'_',1) = species.LOCUSTAG.",
    "cazy_overview": "dbCAN CAZyme calls per protein. species_prefix = species.LOCUSTAG.",
    "gene_info": "Gene coordinates. LOCUSTAG joins species.",
}


def table_name(p: Path) -> str:
    stem = p.name.removesuffix(".parquet")
    return stem.replace(".", "_").replace("-", "_")


def discover(tables: Path, only, include_swissprot):
    exclude = [e for e in DEFAULT_EXCLUDE if not (include_swissprot and e == "swissprot_hits")]
    out = []
    for p in sorted(tables.iterdir()):
        is_file = p.is_file() and p.name.endswith(".parquet")
        is_dir = p.is_dir() and any(p.glob("*.parquet"))
        if not (is_file or is_dir):
            continue
        name = table_name(p)
        if only:
            if name not in only and p.name.removesuffix(".parquet") not in only:
                continue
        elif any(fnmatch.fnmatch(p.name, pat) for pat in exclude):
            continue
        src = str(p) if is_file else str(p / "*.parquet")
        size = p.stat().st_size if is_file else sum(f.stat().st_size for f in p.glob("*.parquet"))
        out.append((name, src, size, p))
    return out


def load_token(token_file: Path):
    if os.environ.get("motherduck_token"):
        return
    if not token_file.exists():
        sys.exit(f"No $motherduck_token and no token file at {token_file}")
    mode = token_file.stat().st_mode
    if mode & (stat.S_IRWXG | stat.S_IRWXO):
        sys.exit(f"{token_file} is readable by group/others; run: chmod 600 {token_file}")
    os.environ["motherduck_token"] = token_file.read_text().strip()


def comment_for(tables: Path, name: str, stem: str):
    """Default description plus build provenance from _manifest/<stem>.json if present."""
    text = DEFAULT_COMMENTS.get(name, "")
    m = tables / "_manifest" / f"{stem}.json"
    if m.exists():
        try:
            d = json.loads(m.read_text())
            text = d.get("description") or text
            prov = ", ".join(f"{k}={d[k]}" for k in ("built_at", "built_by", "schema_version") if k in d)
            if prov:
                text = f"{text} [{prov}]".strip()
        except json.JSONDecodeError:
            pass
    return text or None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tables", type=Path, required=True)
    ap.add_argument("--db", default="bfd")
    ap.add_argument("--only", nargs="*", default=None, help="table names to load (overrides excludes)")
    ap.add_argument("--include-swissprot", action="store_true")
    ap.add_argument("--token-file", type=Path, default=Path("~/.config/motherduck/token").expanduser())
    ap.add_argument("--threads", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 8)))
    ap.add_argument("--memory", default="24GB")
    ap.add_argument("--tmp", default=os.environ.get("SCRATCH", "/tmp") + "/md_duckdb_tmp")
    ap.add_argument("--dry-run", action="store_true", help="list what would load and the local row counts")
    a = ap.parse_args()

    todo = discover(a.tables, a.only, a.include_swissprot)
    total = sum(s for _, _, s, _ in todo)
    print(f"{len(todo)} tables, {total/1e9:.1f} GB parquet -> md:{a.db}", flush=True)

    local = duckdb.connect()
    if a.dry_run:
        for name, src, size, _ in todo:
            n = local.sql(f"SELECT count(*) FROM read_parquet('{src}')").fetchone()[0]
            print(f"  {name:36s} {size/1e9:7.2f} GB {n:>14,d} rows")
        return

    load_token(a.token_file)
    c = duckdb.connect("md:")
    c.sql(f"SET threads={a.threads}; SET memory_limit='{a.memory}'; SET temp_directory='{a.tmp}'")
    c.sql(f"CREATE DATABASE IF NOT EXISTS {a.db}")
    c.sql(f"USE {a.db}")
    existing = {r[0] for r in c.sql(
        f"SELECT table_name FROM duckdb_tables() WHERE database_name='{a.db}'").fetchall()}

    failed = []
    for name, src, size, p in todo:
        want = local.sql(f"SELECT count(*) FROM read_parquet('{src}')").fetchone()[0]
        if name in existing:
            have = c.sql(f"SELECT count(*) FROM {a.db}.main.\"{name}\"").fetchone()[0]
            if have == want:
                print(f"skip  {name:36s} {want:>14,d} rows (already loaded)", flush=True)
                continue
            print(f"redo  {name:36s} have {have:,d} != {want:,d}", flush=True)
        t0 = time.time()
        try:
            c.sql(f"CREATE OR REPLACE TABLE \"{name}\" AS SELECT * FROM read_parquet('{src}')")
            have = c.sql(f"SELECT count(*) FROM \"{name}\"").fetchone()[0]
            if have != want:
                raise RuntimeError(f"row count {have:,d} != source {want:,d}")
            txt = comment_for(a.tables, name, p.name.removesuffix(".parquet"))
            if txt:
                c.sql(f"COMMENT ON TABLE \"{name}\" IS '" + txt.replace("'", "''") + "'")
            dt = time.time() - t0
            print(f"load  {name:36s} {want:>14,d} rows {size/1e9:7.2f} GB {dt:7.0f}s "
                  f"({size/1e6/max(dt,1):.0f} MB/s)", flush=True)
        except Exception as e:  # keep going; resubmit retries only what's missing
            failed.append(name)
            print(f"FAIL  {name}: {str(e)[:300]}", flush=True)

    if failed:
        sys.exit(f"{len(failed)} table(s) failed: {' '.join(failed)} -- resubmit to retry")
    print("done", flush=True)


if __name__ == "__main__":
    main()
