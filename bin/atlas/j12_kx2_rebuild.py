#!/usr/bin/env python3
"""J12: rebuild Analysis E (Kickxellomycota reduction, hypotheses H1-H4) end to end on the local BFD DuckDB.

Runs, in order, the saved SQL of docs/analyses/kx2_kickx_reduction/:
    sql_upstream/  (kw_rep_genomes, kw_pfam_matrix, kx_genomes: inputs that earlier analyses built by hand)
    sql_h1/ ... sql_h4/
into an EMPTY scratch database, and writes every result query and a fingerprint of every kx2_* table to <out>.

Why a driver and not one SQL file: the saved files are the record of what was run in MotherDuck, one statement block at a time. They are
executed here unchanged except for two documented rewrites (SHIMS below): the local views name some columns differently from MotherDuck.

Usage (HPCC, inside a SLURM job; see workflows/atlas/J12_kx2_rebuild.sh):
    j12_kx2_rebuild.py --sql-root docs/analyses/kx2_kickx_reduction --out OUT --scratch-db $SCRATCH/j12.duckdb --tmp $SCRATCH/duckdb_tmp
Options: --stages upstream,h1,h2,h3,h4 (default all); --bfd-db PATH; --threads N; --memory 100GB.
Outputs in OUT:
    run_log.tsv                 stage, file, statement index, seconds, rows returned
    results/<stage>__<file>__q<k>.tsv   one file per statement that returns rows (header + rows)
    table_fingerprint.tsv       per kx2_*/kw_*/kx_* table and column: non-null count and sum (numeric columns), for comparison with MotherDuck
    table_counts.tsv            row count per table
No network, no MotherDuck.
"""
import argparse
import csv
import os
import re
import sys
import time
from pathlib import Path

STAGES = ['upstream', 'h1', 'h2', 'h3', 'h4']
STAGE_DIR = {'upstream': 'sql_upstream', 'h1': 'sql_h1', 'h2': 'sql_h2', 'h3': 'sql_h3', 'h4': 'sql_h4'}

# MotherDuck `bfd.species` has column ORDER_ in the analyses' SQL; the local view is "ORDER". The local asm_stats names two columns differently.
# The shim views live in the scratch database; the SQL text is rewritten to read them. Nothing is written into `bfd`.
SHIM_VIEWS = {
    'j12_species': 'SELECT * EXCLUDE ("ORDER"), "ORDER" AS ORDER_ FROM bfd.species',
    'j12_asm_stats': 'SELECT *, TOTAL_LENGTH AS total_length_bp, GC_PERCENT AS gc_pct FROM bfd.asm_stats',
}
REWRITES = [
    (re.compile(r'\bbfd\.species\b'), 'bfd_scratch.j12_species'),
    (re.compile(r'\bbfd\.asm_stats\b'), 'bfd_scratch.j12_asm_stats'),
]


def rewrite(sql):
    for pat, rep in REWRITES:
        sql = pat.sub(rep, sql)
    return sql


def sql_files(root, stage):
    d = Path(root) / STAGE_DIR[stage]
    return sorted(p for p in d.glob('*.sql'))


def split_statements(con, text):
    """Statements of a SQL file using DuckDB's own parser (safe for semicolons in comments and strings)."""
    return [s.query.strip() for s in con.extract_statements(text) if s.query.strip()]


def write_rows(path, cols, rows):
    with open(path, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(cols)
        for r in rows:
            w.writerow(['' if v is None else v for v in r])


def run_stage(con, root, stage, out, log):
    for p in sql_files(root, stage):
        text = rewrite(p.read_text())
        for k, stmt in enumerate(split_statements(con, text), 1):
            t0 = time.time()
            cur = con.execute(stmt)
            nrows = ''
            if cur.description is not None:
                cols = [d[0] for d in cur.description]
                rows = cur.fetchall()
                nrows = len(rows)
                if not (len(cols) == 1 and cols[0] == 'Count' and nrows == 1):   # skip the "N rows affected" echo of DDL/DML
                    write_rows(out / 'results' / f'{stage}__{p.stem}__q{k}.tsv', cols, rows)
            dt = time.time() - t0
            log.writerow([stage, p.name, k, f'{dt:.1f}', nrows])
            sys.stdout.flush()
            print(f'{stage} {p.name} stmt {k}: {dt:.1f}s rows={nrows}', flush=True)


def fingerprint(con, out):
    tabs = [r[0] for r in con.execute(
        "SELECT table_name FROM duckdb_tables() WHERE database_name='bfd_scratch' ORDER BY 1").fetchall()]
    with open(out / 'table_counts.tsv', 'w', newline='') as fc, open(out / 'table_fingerprint.tsv', 'w', newline='') as ff:
        wc = csv.writer(fc, delimiter='\t', lineterminator='\n')
        wf = csv.writer(ff, delimiter='\t', lineterminator='\n')
        wc.writerow(['table', 'n_rows'])
        wf.writerow(['table', 'column', 'type', 'n_nonnull', 'sum'])
        for t in tabs:
            n = con.execute(f'SELECT count(*) FROM bfd_scratch."{t}"').fetchone()[0]
            wc.writerow([t, n])
            for col, typ in con.execute(
                    "SELECT column_name, data_type FROM duckdb_columns() WHERE database_name='bfd_scratch' AND table_name=? ORDER BY column_index", [t]).fetchall():
                if re.match(r'(TINYINT|SMALLINT|INTEGER|BIGINT|HUGEINT|UBIGINT|UINTEGER|FLOAT|DOUBLE|DECIMAL)', typ):
                    nn, sm = con.execute(f'SELECT count("{col}"), sum("{col}"::DOUBLE) FROM bfd_scratch."{t}"').fetchone()
                    wf.writerow([t, col, typ, nn, '' if sm is None else f'{sm:.6g}'])
                else:
                    nn = con.execute(f'SELECT count("{col}") FROM bfd_scratch."{t}"').fetchone()[0]
                    wf.writerow([t, col, typ, nn, ''])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--sql-root', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--scratch-db', required=True)
    ap.add_argument('--bfd-db', default=None)
    ap.add_argument('--lib', default=None, help='directory holding bfd_connect.py (default: <repo>/lib)')
    ap.add_argument('--tmp', default=None)
    ap.add_argument('--threads', type=int, default=16)
    ap.add_argument('--memory', default='100GB')
    ap.add_argument('--stages', default=','.join(STAGES))
    a = ap.parse_args(argv)

    lib = a.lib or str(Path(__file__).resolve().parents[2] / 'lib')
    sys.path.insert(0, lib)
    import bfd_connect

    out = Path(a.out)
    (out / 'results').mkdir(parents=True, exist_ok=True)
    con = bfd_connect.connect(a.bfd_db, a.scratch_db, threads=a.threads, memory=a.memory, tmp=a.tmp)
    for name, body in SHIM_VIEWS.items():
        con.execute(f'CREATE OR REPLACE VIEW bfd_scratch.{name} AS {body}')
    con.execute('USE bfd_scratch')           # bare names in the saved SQL are scratch tables, as in MotherDuck
    stages = [s for s in a.stages.split(',') if s]
    bad = [s for s in stages if s not in STAGES]
    if bad:
        sys.exit(f'unknown stage(s): {bad}')
    with open(out / 'run_log.tsv', 'a', newline='') as lf:
        log = csv.writer(lf, delimiter='\t', lineterminator='\n')
        if lf.tell() == 0:
            log.writerow(['stage', 'file', 'stmt', 'seconds', 'rows'])
        for s in stages:
            run_stage(con, a.sql_root, s, out, log)
    fingerprint(con, out)
    print('J12 done', flush=True)


if __name__ == '__main__':
    main()
