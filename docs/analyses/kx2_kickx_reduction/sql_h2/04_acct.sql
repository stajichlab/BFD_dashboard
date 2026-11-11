-- 04_acct.sql
-- Purpose: per-genome accounting: assembly = sum(gene span) - overlap + sum(positive gaps) + remainder.
--   Gene span already contains CDS + introns (UTR = 0, see 02_tx.sql). remainder = contig-edge flanks, gene-free contigs, N.
--   masked_bases (from 01) overlaps the intergenic/intron terms and is reported as a memo item, not additive.
-- Inputs: bfd_scratch.kx2_h2_genomes, kx2_h2_tx, kx2_h2_ig.
-- Output: bfd_scratch.kx2_h2_acct
-- Report numbers: feeds every table (lineage Mb terms, per-gene intergenic, fraction of genes with introns, remainder).
CREATE OR REPLACE TABLE kx2_h2_acct AS
WITH t AS (SELECT genome, count(*) n_genes, sum(cds) cds_bp, sum(span-exon_len) intron_bp, sum(exon_len-cds) utr_bp, sum(span) span_bp,
  avg((n_exon>1)::int) f_intron, avg(cds) mean_cds, avg(span) mean_span, avg(n_exon-1) ipt, avg(is_partial::int) f_partial FROM kx2_h2_tx GROUP BY 1),
i AS (SELECT genome, count(*) n_gaps, sum(greatest(distance,0)) gap_bp, sum(greatest(-distance,0)) overlap_bp, avg((distance<=0)::int) f_ov, avg((distance<100)::int) f_lt100,
  median(distance) ig_med, quantile_cont(distance,0.25) ig_q25 FROM kx2_h2_ig GROUP BY 1)
SELECT g.*, t.*EXCLUDE(genome), i.*EXCLUDE(genome),
  g.total_length_bp - t.span_bp - i.gap_bp + i.overlap_bp AS remainder_bp
FROM kx2_h2_genomes g JOIN t ON t.genome=g.LOCUSTAG JOIN i ON i.genome=g.LOCUSTAG;
