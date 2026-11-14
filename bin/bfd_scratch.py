#!/usr/bin/env python3
"""Manage the local bfd_scratch database (replacement for MotherDuck bfd_scratch).

    bfd_scratch.py init                 create the scratch file (does not need BFD.duckdb)
    bfd_scratch.py list                 scratch tables by analysis prefix, with row estimates
    bfd_scratch.py sql "SELECT ..."     run SQL with bfd (read-only) and bfd_scratch attached
    bfd_scratch.py sql -f query.sql     same, from a file
    bfd_scratch.py drop kw_             drop every table whose name starts with kw_ (asks first)

In SQL, write bfd_scratch.<prefix>_<name> for scratch tables and bfd.<table> (or the bare
table name) for BFD data. Tables in bfd_scratch need an analysis prefix.

Paths: --bfd-db / $BFD_DB and --scratch-db / $BFD_SCRATCH_DB (see lib/bfd_connect.py).
Run with /usr/bin/python3.12.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "lib"))
import bfd_connect as bc  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bfd-db", type=Path)
    ap.add_argument("--scratch-db", type=Path)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("init")
    sub.add_parser("list")
    q = sub.add_parser("sql")
    q.add_argument("query", nargs="?")
    q.add_argument("-f", "--file", type=Path)
    d = sub.add_parser("drop")
    d.add_argument("prefix")
    d.add_argument("--yes", action="store_true")
    a = ap.parse_args()

    if a.cmd == "init":
        import duckdb
        import os
        sp = Path(a.scratch_db or os.environ.get("BFD_SCRATCH_DB") or bc.DEFAULT_SCRATCH_DB)
        sp.parent.mkdir(parents=True, exist_ok=True)
        duckdb.connect(str(sp)).close()
        print(f"scratch database ready: {sp}")
        return

    con = bc.connect(a.bfd_db, a.scratch_db)

    if a.cmd == "list":
        rows = bc.scratch_tables(con)
        if not rows:
            print("bfd_scratch is empty")
        for name, prefix, est in rows:
            print(f"{prefix:14s} {name:48s} ~{est:,d} rows")
    elif a.cmd == "sql":
        text = a.file.read_text() if a.file else a.query
        if not text:
            sys.exit("give a query or -f file")
        res = con.execute(text)
        if res.description:
            print(res.fetchdf().to_string(index=False))
    elif a.cmd == "drop":
        bc.check_prefix(a.prefix if a.prefix.endswith("_") else a.prefix + "_")
        victims = [n for n, _, _ in bc.scratch_tables(con) if n.startswith(a.prefix)]
        if not victims:
            sys.exit(f"no scratch tables start with {a.prefix}")
        print("will drop:", *victims, sep="\n  ")
        if not a.yes and input("type 'yes' to drop: ").strip() != "yes":
            sys.exit("cancelled")
        for n in victims:
            con.execute(f'DROP TABLE bfd_scratch."{n}"')
        print(f"dropped {len(victims)} table(s)")


if __name__ == "__main__":
    main()
