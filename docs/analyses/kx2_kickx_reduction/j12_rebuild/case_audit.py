"""Audit of case-sensitivity effects in the shared BFD parquet tables (read-only).
Usage: case_audit.py <tables_dir>
Prints: schema of the gene_* tables; per-genome partial-transcript fraction histogram; lower-case counts in sequence-derived text columns.
"""
import sys
import duckdb

T = sys.argv[1]
c = duckdb.connect()
c.execute("SET threads=8")
for t in ['gene_transcripts', 'gene_introns', 'gene_proteins', 'gene_CDS', 'gene_exons', 'gene_info', 'gene_trnas']:
    cols = c.execute(f"DESCRIBE SELECT * FROM read_parquet('{T}/{t}.parquet')").fetchall()
    print(t, [(r[0], r[1]) for r in cols])
print()
print('== per-genome partial fraction (gene_transcripts) ==')
q = f"""
WITH g AS (SELECT split_part(transcript_id,'_',1) lt, count(*) n, avg(is_partial::int) fp,
                  avg((has_start_codon='FALSE')::int) fstart, avg((has_stop_codon='FALSE')::int) fstop
           FROM read_parquet('{T}/gene_transcripts.parquet') GROUP BY 1)
SELECT count(*) genomes, count(*) FILTER (WHERE fp>=0.5) ge50, count(*) FILTER (WHERE fp>=0.2 AND fp<0.5) b20_50,
       count(*) FILTER (WHERE fp>=0.1 AND fp<0.2) b10_20, count(*) FILTER (WHERE fp>=0.05 AND fp<0.1) b5_10,
       count(*) FILTER (WHERE fp<0.05) lt5, round(median(fp),4) median_fp, round(max(fp),3) max_fp, sum(n) tx FROM g
"""
print(c.execute(q).fetchall())
print('genomes with partial fraction >= 0.10 (first 25):')
print(c.execute(f"""
SELECT split_part(transcript_id,'_',1) lt, count(*) n, round(avg(is_partial::int),3) fp,
       round(avg((has_start_codon='FALSE')::int),3) fstart, round(avg((has_stop_codon='FALSE')::int),3) fstop
FROM read_parquet('{T}/gene_transcripts.parquet') GROUP BY 1 HAVING avg(is_partial::int) >= 0.10 ORDER BY fp DESC LIMIT 25""").fetchall())
print()
print('== lower-case in sequence-derived columns ==')
print('gene_transcripts distinct has_start_codon / has_stop_codon values:',
      c.execute(f"SELECT has_start_codon, has_stop_codon, count(*) FROM read_parquet('{T}/gene_transcripts.parquet') GROUP BY 1,2 ORDER BY 3 DESC").fetchall())
try:
    print('gene_introns splice_5 lower-case rows:', c.execute(
        f"SELECT count(*) FILTER (WHERE splice_5 <> upper(splice_5)), count(*) FROM read_parquet('{T}/gene_introns.parquet')").fetchall())
    print('gene_introns splice_3 lower-case rows:', c.execute(
        f"SELECT count(*) FILTER (WHERE splice_3 <> upper(splice_3)), count(*) FROM read_parquet('{T}/gene_introns.parquet')").fetchall())
    print('gene_introns top splice_5 / splice_3 pairs:', c.execute(
        f"SELECT splice_5, splice_3, count(*) FROM read_parquet('{T}/gene_introns.parquet') GROUP BY 1,2 ORDER BY 3 DESC LIMIT 8").fetchall())
except Exception as e:
    print('gene_introns check failed:', e)
