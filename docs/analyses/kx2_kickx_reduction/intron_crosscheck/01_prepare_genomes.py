#!/usr/bin/env python3
"""Add assembly ID, contig count, N50 and BUSCO to genomes.tsv from the BFD parquet tables and check the assembly file exists.

Usage (HPCC, python with duckdb, e.g. /usr/bin/python3.12): 01_prepare_genomes.py genomes.tsv genomes_resolved.tsv
Read-only on the BFD tables. Exits non-zero if any genome has no assembly ID or no assembly file.
"""
import csv
import os
import sys

import duckdb

T = '/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables'
CLEAN = '/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes'
src, dst = sys.argv[1], sys.argv[2]
rows = list(csv.DictReader(open(src), delimiter='\t'))
con = duckdb.connect()
res = con.execute(f"""
  SELECT s.LOCUSTAG, s.ASMID, a.contig_count, a.N50_bp, b.complete_pct, b.duplicated_pct
  FROM '{T}/species.parquet' s JOIN '{T}/asm_stats.parquet' a USING (ASMID) JOIN '{T}/busco_genome.parquet' b USING (ASMID)
""").fetchall()
info = {r[0]: r for r in res}
bad = []
with open(dst, 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t', lineterminator='\n')
    w.writerow(['locustag', 'genus', 'ord', 'phylum', 'set', 'role', 'asmid', 'contig_count', 'n50', 'busco_c', 'busco_d', 'assembly'])
    for r in rows:
        i = info.get(r['locustag'])
        if i is None:
            bad.append((r['locustag'], 'no species/asm_stats/busco row'))
            continue
        asm = f"{CLEAN}/{i[1]}.fa.gz"
        if not (os.path.exists(asm) and os.path.getsize(asm) > 0):
            bad.append((r['locustag'], f'missing or empty {asm}'))
        w.writerow([r['locustag'], r['genus'], r['ord'], r['phylum'], r['set'], r['role'], i[1], i[2], i[3], i[4], i[5], asm])
print(f'{len(rows)} genomes written to {dst}; problems: {bad}')
sys.exit(1 if bad else 0)
