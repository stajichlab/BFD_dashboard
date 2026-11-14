-- 03_swissprot_best.sql
-- Purpose : SwissProt best-hit (name, organism, identity, coverage, e-value) for every protein of the 169 genomes.
-- Inputs  : bfd.swissprot_best (identity/coverage stored x10), bfd.swissprot_subject (gene/protein names), kx2_ku_genomes
-- Output  : bfd_scratch.kx2_ku_sp  (873,235 rows when first run)
-- Thresholds: none here (name rules in 05 apply pident>=25, scov>=40, -log10 e>=10).
-- Reports : none directly (intermediate).
CREATE OR REPLACE TABLE bfd_scratch.kx2_ku_sp AS
SELECT b.protein_id, split_part(b.protein_id,'_',1) AS locustag, b.accession, s.gene_name, s.protein_name, s.organism,
       b.pident_x10/10.0 AS pident, b.qcov_x10/10.0 AS qcov, b.scov_x10/10.0 AS scov, b.bitscore, b.neg_log10_evalue AS nlog10e, b.q_len, b.s_len
FROM bfd.swissprot_best b JOIN bfd.swissprot_subject s ON s.sp_id=b.sp_id
WHERE split_part(b.protein_id,'_',1) IN (SELECT locustag FROM bfd_scratch.kx2_ku_genomes);
