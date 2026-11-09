#!/usr/bin/env python3
"""Run one atlas SQL job file with DuckDB. Restartable: every output block is skipped when
its file exists, written to <name>.tmp and renamed when done.

SQL file format:
  statements before the first '-- @output' marker are a preamble (views, settings) and always run;
  '-- @output <relative/path.parquet>' starts a block; the block's statements run with the
  placeholder {OUT} replaced by the temporary output path.
Placeholders available everywhere: {T} tables dir, {A} atlas_data dir, {OUT} (blocks only).
"""
import argparse
import hashlib
import json
import os
import re
import sys
import time

import duckdb


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--sql', required=True)
    ap.add_argument('--tables', default='/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables')
    ap.add_argument('--atlas', default='/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/atlas_data')
    ap.add_argument('--tmp', required=True, help='DuckDB temp dir (node-local $SCRATCH)')
    ap.add_argument('--threads', type=int, default=16)
    ap.add_argument('--memory', default='100GB')
    ap.add_argument('--force', action='store_true')
    a = ap.parse_args()

    os.makedirs(a.tmp, exist_ok=True)
    text = open(a.sql).read()
    sha = hashlib.sha256(text.encode()).hexdigest()[:12]
    parts = re.split(r'^-- @output\s+(\S+)\s*$', text, flags=re.M)
    preamble = parts[0]
    blocks = list(zip(parts[1::2], parts[2::2]))

    con = duckdb.connect()
    con.execute(f"SET threads={a.threads}; SET memory_limit='{a.memory}'; SET temp_directory='{a.tmp}'; "
                "SET preserve_insertion_order=false; SET enable_progress_bar=false;")
    sub = lambda s, out='': s.replace('{T}', a.tables).replace('{A}', a.atlas).replace('{OUT}', out)
    t0 = time.time()
    print(f'[{time.strftime("%H:%M:%S")}] preamble', flush=True)
    con.execute(sub(preamble))
    for name, body in blocks:
        final = os.path.join(a.atlas, name)
        if os.path.exists(final) and not a.force:
            print(f'[{time.strftime("%H:%M:%S")}] skip {name} (exists)', flush=True)
            continue
        os.makedirs(os.path.dirname(final), exist_ok=True)
        tmp = final + '.tmp'
        if os.path.exists(tmp):
            os.remove(tmp)
        t1 = time.time()
        print(f'[{time.strftime("%H:%M:%S")}] start {name}', flush=True)
        try:
            con.execute(sub(body, tmp))
        except Exception as e:
            print(f'FAILED {name}: {e}', file=sys.stderr, flush=True)
            if os.path.exists(tmp):
                os.remove(tmp)
            raise
        os.replace(tmp, final)
        n = con.execute(f"SELECT count(*) FROM '{final}'").fetchone()[0]
        man = {'output': name, 'rows': n, 'sql': os.path.basename(a.sql), 'sql_sha256_12': sha,
               'built_at': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), 'seconds': round(time.time() - t1, 1),
               'slurm_job': os.environ.get('SLURM_JOB_ID')}
        with open(final + '.manifest.json', 'w') as o:
            json.dump(man, o, indent=1)
        print(f'[{time.strftime("%H:%M:%S")}] done {name} rows={n} {man["seconds"]}s', flush=True)
    print(f'all done in {time.time() - t0:.0f}s', flush=True)


if __name__ == '__main__':
    main()
