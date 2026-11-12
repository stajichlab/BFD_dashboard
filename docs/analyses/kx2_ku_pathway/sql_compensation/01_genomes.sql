-- Purpose: genome table (taxonomy, BUSCO, assembly stats, gene counts) for the Ku-loss compensation analysis (Q2).
-- Run with default database bfd_scratch; reads bfd.* (read-only).
-- Inputs: bfd.species, bfd.busco_genome, bfd.asm_stats, bfd.gene_proteins
-- Output: bfd_scratch.kx2_kc_genomes (one row per LOCUSTAG, 23,680 genomes)
-- Thresholds: none here (quality flags are set in 06_wide.sql: q_loose, q_strict).
-- Produces: genome universe and n_genes (proteome-size control) used everywhere.
CREATE OR REPLACE TABLE kx2_kc_genomes AS
WITH np AS (SELECT split_part(protein_id,'_',1) AS LOCUSTAG, count(distinct gene_id) n_genes, count(*) n_prot FROM bfd.gene_proteins GROUP BY 1)
SELECT s.LOCUSTAG, s.ASMID, s.PHYLUM, s.SUBPHYLUM, s.CLASS, s."ORDER" AS ORD, s.FAMILY, s.GENUS, s.SPECIES, s.STRAIN,
 b.complete_pct, b.duplicated_pct, a.N50_bp, a.contig_count, a.total_length_bp, a.t2t_scaffolds, np.n_genes, np.n_prot
FROM bfd.species s LEFT JOIN bfd.busco_genome b USING (ASMID) LEFT JOIN bfd.asm_stats a USING (ASMID) LEFT JOIN np ON np.LOCUSTAG=s.LOCUSTAG;
-- check: 23680 rows, 23451 with n_genes, 23582 with BUSCO/N50
SELECT count(*), count(n_genes), count(complete_pct), count(N50_bp) FROM kx2_kc_genomes;
