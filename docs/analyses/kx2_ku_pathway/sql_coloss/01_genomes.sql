-- 01_genomes.sql
-- Purpose : genome set for the Ku-pathway survey = 63 backbone genus representatives
--           + all Coemansia strains passing QC, with QC covariates and group labels.
-- Inputs  : bfd_scratch.kx2_tree_taxa (tree='backbone'), bfd.species, bfd.busco_genome,
--           bfd.asm_stats, bfd.gene_proteins
-- Output  : bfd_scratch.kx2_ku_genomes (one row per genome; set_name = 'backbone' | 'coemansia')
-- Thresholds: Coemansia strains: BUSCO complete >= 60 AND n_proteins >= 3000.
--           Backbone taxa kept regardless of QC (flagged via covariates; caveat taxa Tieghemiomyces, Spiromyces).
-- Reports : n genomes per group; QC covariates (busco_complete, busco_dup, n50, n_proteins).
-- Note    : the backbone Coemansia representative F0DE0BF4 is also in the Coemansia set if it passes QC;
--           set_name is then 'backbone' (row kept once), flag in_coemansia_set marks it.
CREATE OR REPLACE TABLE bfd_scratch.kx2_ku_genomes AS
WITH np AS (
  SELECT split_part(protein_id,'_',1) AS locustag, count(*) AS n_proteins
  FROM bfd.gene_proteins GROUP BY 1
),
cand AS (
  SELECT t.locustag, 'backbone' AS set_name, t.role, t.genus AS genus_t, t.ord AS ord_t, t.phylum AS phylum_t
  FROM bfd_scratch.kx2_tree_taxa t WHERE t.tree='backbone'
  UNION ALL
  SELECT s.LOCUSTAG, 'coemansia', 'ingroup', 'Coemansia', 'Kickxellales', 'Kickxellomycota'
  FROM bfd.species s JOIN bfd.busco_genome b USING(ASMID) JOIN np ON np.locustag=s.LOCUSTAG
  WHERE s.GENUS='Coemansia' AND b.complete_pct>=60 AND np.n_proteins>=3000
),
dedup AS (
  SELECT locustag, arg_min(set_name, CASE set_name WHEN 'backbone' THEN 0 ELSE 1 END) AS set_name,
         arg_min(role, CASE set_name WHEN 'backbone' THEN 0 ELSE 1 END) AS role,
         first(genus_t) AS genus, first(ord_t) AS ord, first(phylum_t) AS phylum,
         bool_or(set_name='coemansia') AS in_coemansia_set
  FROM cand GROUP BY locustag
)
SELECT d.locustag, s.ASMID, s.SPECIES, s.STRAIN, d.set_name, d.role, d.genus, d.ord, d.phylum, d.in_coemansia_set,
  CASE
    WHEN d.genus IN ('Coemansia','Kickxella','Linderina','Dipsacomyces','Martensiomyces') THEN 'Kickxellales_core'
    WHEN d.genus IN ('Mycoemilia','Spiromyces') THEN 'Mycoemilia_Spiromyces'
    WHEN d.ord='Harpellales' THEN 'Harpellales'
    WHEN d.ord='Dimargaritales' THEN 'Dimargaritales'
    WHEN d.genus='Ramicandelaber' THEN 'Ramicandelaber'
    WHEN d.ord='Zoopagales' THEN 'Zoopagales'
    WHEN d.ord='Entomophthorales' OR d.genus='Basidiobolus' THEN 'Entomophthorales_Basidiobolus'
    WHEN d.role='mucoromycot' THEN 'Mucoromycota_s.l.'
    WHEN d.role='root' THEN 'Root'
    ELSE 'other' END AS grp,
  b.complete_pct AS busco_complete, b.duplicated_pct AS busco_dup, b.missing_pct AS busco_missing,
  a.N50_bp AS n50, a.total_length_bp, np.n_proteins
FROM dedup d
JOIN bfd.species s ON s.LOCUSTAG=d.locustag
LEFT JOIN bfd.busco_genome b ON b.ASMID=s.ASMID
LEFT JOIN bfd.asm_stats a ON a.ASMID=s.ASMID
LEFT JOIN np ON np.locustag=d.locustag;
