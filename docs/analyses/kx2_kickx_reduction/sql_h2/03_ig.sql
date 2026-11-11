-- 03_ig.sql
-- Purpose: adjacent-gene intergenic distances annotated with orientation.
-- Inputs: bfd_scratch.kx2_h2_genomes, bfd.gene_info (strand), bfd.gene_intergenic_distances (distance = right.start - left.end).
-- Output: bfd_scratch.kx2_h2_ig
-- Orientation: left +, right - = convergent; left -, right + = divergent; otherwise tandem.
-- Note: distance is >= 1 everywhere except a few tandem pairs (min -6): gene overlap is effectively absent from the table.
-- Report numbers: orientation-specific intergenic results; overlap fraction (0).
CREATE OR REPLACE TABLE kx2_h2_ig AS
WITH gl AS (SELECT LOCUSTAG FROM kx2_h2_genomes WHERE lineage IS NOT NULL),
s AS (SELECT gene_id, strand FROM bfd.gene_info WHERE LOCUSTAG IN (SELECT LOCUSTAG FROM gl))
SELECT d.species_prefix AS genome, d.distance,
 CASE WHEN l.strand=1 AND r.strand=-1 THEN 'convergent' WHEN l.strand=-1 AND r.strand=1 THEN 'divergent' ELSE 'tandem' END AS orient
FROM bfd.gene_intergenic_distances d JOIN s l ON l.gene_id=d.left_gene JOIN s r ON r.gene_id=d.right_gene
WHERE d.species_prefix IN (SELECT LOCUSTAG FROM gl);
