-- 09_keydomain_sensitivity.sql
-- Purpose : sensitivity analysis. The multi-domain rules in 04 can under-call (a protein can lose one of two
--           required domains to annotation noise). Re-score each Pfam-defined gene with ONE key domain
--           (any of the listed accessions) at i-evalue <= 1e-3, and compare with the strict rule.
-- Inputs  : bfd_scratch.kx2_ku_pfam (02), kx2_ku_genomes, kx2_ku_matrix
-- Output  : bfd_scratch.kx2_ku_key (gene, acc), bfd_scratch.kx2_ku_keymatrix (gene x genome, key_present),
--           bfd_scratch.kx2_ku_key_prev (gene x group, strict vs key-domain counts)
-- Thresholds: i-evalue <= 1e-3 (kx2_ku_pfam keeps <= 1e-3), any one key accession.
-- Reports : relaxed presence by group; relaxed false-absence rate of controls; confirms whether absences
--           (Ku, DNA-PKcs, XLF, Pol4 ...) survive the relaxed rule.
-- Note    : gene names must match kx2_ku_rules.gene.
CREATE OR REPLACE TABLE bfd_scratch.kx2_ku_key AS
SELECT * FROM (VALUES
 ('Ku','PF02735'),('Ku','PF03731'),('Ku','PF03730'),
 ('DNA-PKcs','PF08163'),('DNA-PKcs','PF20500'),('DNA-PKcs','PF20502'),('DNA-PKcs','PF19704'),
 ('Lig4','PF11411'),
 ('XRCC4/Lif1','PF06632'),('XLF/Nej1','PF09302'),('PAXX','PF15384'),
 ('Pol4 (X-family Pol mu/lambda)','PF14792'),('Aprataxin (Hnt3)','PF16278'),('TDP1','PF06087'),
 ('Mre11','PF04152'),('Rad50','PF04423'),('Sae2/CtIP','PF08573'),('Dna2','PF08696'),
 ('ATM/Tel1','PF11640'),('ATR/Mec1','PF08064'),('TERT','PF12009'),('Est1','PF10374'),('Rif1','PF12231'),
 ('Taz1/TRF','PF08558'),('Stn1','PF10451'),('Ten1','PF15490'),
 ('CTRL Rad51 family','PF08423'),('CTRL Pol delta (Pol3)','PF24055'),('CTRL Pol epsilon (Pol2)','PF08490'),('CTRL Mcm2','PF12619')
) t(gene, acc);

CREATE OR REPLACE TABLE bfd_scratch.kx2_ku_keymatrix AS
SELECT k.gene, g.locustag, g.genus, g.grp, g.set_name,
       EXISTS (SELECT 1 FROM bfd_scratch.kx2_ku_pfam f WHERE f.locustag=g.locustag AND f.acc IN (SELECT acc FROM bfd_scratch.kx2_ku_key k2 WHERE k2.gene=k.gene)) AS key_present
FROM (SELECT DISTINCT gene FROM bfd_scratch.kx2_ku_key) k CROSS JOIN bfd_scratch.kx2_ku_genomes g;

CREATE OR REPLACE TABLE bfd_scratch.kx2_ku_key_prev AS
WITH u AS (
  SELECT gene, locustag, key_present, grp AS grp2 FROM bfd_scratch.kx2_ku_keymatrix WHERE set_name='backbone'
  UNION ALL
  SELECT gene, locustag, key_present, 'Coemansia_strains' FROM bfd_scratch.kx2_ku_keymatrix WHERE genus='Coemansia')
SELECT gene, grp2, count(*) AS n_genomes, sum(key_present::INT) AS n_key_present FROM u GROUP BY ALL;
