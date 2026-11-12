-- 02_tx.sql
-- Purpose: per-gene longest transcript with span, CDS length and exon length (for intron and UTR accounting).
-- Inputs: bfd_scratch.kx2_h2_genomes (lineage NOT NULL), bfd.gene_transcripts, bfd.gene_exons.
-- Output: bfd_scratch.kx2_h2_tx (12,423,748 rows)
-- Logic: one transcript per gene (max CDS_length, ties by transcript_id). intron_bp = span - exon_len; UTR bp = exon_len - cds.
--   Genome key = split_part(gene_id,'_',1) = LOCUSTAG. Big tables are filtered to chosen genomes first.
-- Report numbers: gene counts, CDS/intron bp, fraction of genes with introns, and the UTR check
--   (exon_len - cds > 3 is never true in any lineage: UTR is not annotated, so the UTR term is 0 by construction).
CREATE OR REPLACE TABLE kx2_h2_tx AS
WITH gl AS (SELECT LOCUSTAG FROM kx2_h2_genomes WHERE lineage IS NOT NULL),
t AS (SELECT split_part(gene_id,'_',1) AS genome, gene_id, transcript_id, start, "end", CDS_length, is_partial,
   row_number() OVER (PARTITION BY gene_id ORDER BY CDS_length DESC, transcript_id) AS rn
   FROM bfd.gene_transcripts WHERE split_part(gene_id,'_',1) IN (SELECT LOCUSTAG FROM gl)),
t1 AS (SELECT * FROM t WHERE rn=1),
e AS (SELECT transcript_id, sum("end"-start+1) AS exon_len, count(*) AS n_exon FROM bfd.gene_exons
   WHERE split_part(transcript_id,'_',1) IN (SELECT LOCUSTAG FROM gl) GROUP BY 1)
SELECT t1.genome, t1.gene_id, t1.transcript_id, t1."end"-t1.start+1 AS span, t1.CDS_length AS cds, e.exon_len, e.n_exon, t1.is_partial
FROM t1 LEFT JOIN e USING(transcript_id);

-- UTR annotation check (run afterwards; all rows returned 0):
-- SELECT lineage, count(*) ng, round(median(f_utr),3) med_frac_tx_utr, round(quantile_cont(f_utr,0.9),3) q90, round(median(mean_utr),1) med_mean_utr, sum((f_utr>0.2)::int) n_genomes_utr_annot
-- FROM (SELECT g.lineage, g.LOCUSTAG, avg((t.exon_len-t.cds>3)::int) f_utr, avg(t.exon_len-t.cds) mean_utr FROM kx2_h2_tx t JOIN kx2_h2_genomes g ON g.LOCUSTAG=t.genome GROUP BY 1,2) GROUP BY 1 ORDER BY 1;
