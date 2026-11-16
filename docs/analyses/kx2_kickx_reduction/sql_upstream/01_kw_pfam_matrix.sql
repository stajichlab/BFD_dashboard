-- 01_kw_pfam_matrix.sql (J12 rebuild of an input table)
-- Purpose: genome x Pfam matrix over the species representatives. Gate domain_i_evalue <= 1e-5, version suffix stripped.
-- Origin: kingdom_domains/README.md (same statement).
-- Output: bfd_scratch.kw_pfam_matrix  (MotherDuck: 32,265,977 rows, 7,392 genomes x 15,052 Pfams)
CREATE OR REPLACE TABLE bfd_scratch.kw_pfam_matrix AS
SELECT split_part(protein_id,'_',1) AS genome, split_part(pfam_acc,'.',1) AS pfam_acc, any_value(pfam_id) AS pfam_id,
       count(DISTINCT protein_id) AS n_proteins, count(*) AS n_domains
FROM bfd.pfam
WHERE domain_i_evalue <= 1e-5 AND split_part(protein_id,'_',1) IN (SELECT LOCUSTAG FROM bfd_scratch.kw_rep_genomes)
GROUP BY 1, 2;
-- check
SELECT count(*) AS n_rows, count(DISTINCT genome) AS n_genomes, count(DISTINCT pfam_acc) AS n_pfams FROM bfd_scratch.kw_pfam_matrix;
