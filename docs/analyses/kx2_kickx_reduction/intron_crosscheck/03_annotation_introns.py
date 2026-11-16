#!/usr/bin/env python3
"""Introns per transcript from the funannotate models in BFD, for the genomes in genomes_resolved.tsv.

Usage (HPCC, /usr/bin/python3.12, via 03_annotation_introns.sh): 03_annotation_introns.py genomes_resolved.tsv annotation_introns.tsv
Reads gene_transcripts.parquet and gene_introns.parquet (read-only). Genome = split_part(transcript_id, '_', 1).
Output columns: locustag, n_transcripts, n_introns, introns_per_transcript, frac_transcripts_with_intron.
Note: this is a genome-wide value over all annotated transcripts; the miniprot metric covers 1,122 conserved single-copy genes
only, so compare ranks and relative group order, not absolute values.
"""
import csv
import sys

import duckdb

T = '/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables'
rows = list(csv.DictReader(open(sys.argv[1]), delimiter='\t'))
locs = [r['locustag'] for r in rows]
con = duckdb.connect()
con.execute("SET threads=4; SET memory_limit='24GB'")
con.execute("CREATE TABLE g(loc VARCHAR)")
con.executemany("INSERT INTO g VALUES (?)", [(x,) for x in locs])
res = con.execute(f"""
  WITH tx AS (SELECT split_part(transcript_id,'_',1) AS loc, transcript_id FROM '{T}/gene_transcripts.parquet'
              WHERE split_part(transcript_id,'_',1) IN (SELECT loc FROM g)),
       ic AS (SELECT transcript_id, count(*) AS n FROM '{T}/gene_introns.parquet'
              WHERE split_part(transcript_id,'_',1) IN (SELECT loc FROM g) GROUP BY 1)
  SELECT tx.loc, count(*) AS n_tx, sum(coalesce(ic.n,0)) AS n_introns, count(*) FILTER (WHERE ic.n > 0) AS n_tx_with_intron
  FROM tx LEFT JOIN ic USING (transcript_id) GROUP BY 1 ORDER BY 1
""").fetchall()
with open(sys.argv[2], 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t', lineterminator='\n')
    w.writerow(['locustag', 'n_transcripts', 'n_introns', 'introns_per_transcript', 'frac_transcripts_with_intron'])
    for loc, n_tx, n_i, n_w in res:
        w.writerow([loc, n_tx, n_i, f'{n_i / n_tx:.4f}', f'{n_w / n_tx:.4f}'])
missing = sorted(set(locs) - {r[0] for r in res})
print(f'{len(res)} genomes with transcripts; missing: {missing}')
sys.exit(1 if missing else 0)
