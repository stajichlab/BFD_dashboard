-- 08_prevalence_controls.sql
-- Purpose : (a) presence prevalence per gene per group; (b) false-absence estimate from ubiquitous controls.
-- Inputs  : bfd_scratch.kx2_ku_matrix, kx2_ku_cov
-- Output  : bfd_scratch.kx2_ku_prev      (gene, conf, category, grp2, n_genomes, n_present, frac)
--           bfd_scratch.kx2_ku_ctrl_cov  (per genome: n strict controls found of 4, QC covariates)
--           bfd_scratch.kx2_ku_ctrl_grp  (per group: mean control detection = 1 - false-absence rate estimate)
-- Group definitions (grp2): backbone-only groups from kx2_ku_genomes.grp for the 63 backbone genera (one per genus,
--           so genera, not strains, are the unit); 'Coemansia_strains' = all 107 QC-passing Coemansia genomes
--           (including the backbone Coemansia representative).
-- Thresholds: strict controls = rules with is_ctrl=true (Rad51 family, Pol delta, Pol epsilon, Mcm2).
--           False-absence estimate for a group = 1 - mean(fraction of the 4 controls detected).
-- Reports : presence-by-group table; false-absence rate by group.
CREATE OR REPLACE TABLE bfd_scratch.kx2_ku_prev AS
WITH u AS (
  SELECT gene, conf, category, locustag, present, grp AS grp2 FROM bfd_scratch.kx2_ku_matrix WHERE set_name='backbone'
  UNION ALL
  SELECT gene, conf, category, locustag, present, 'Coemansia_strains' FROM bfd_scratch.kx2_ku_matrix WHERE genus='Coemansia'
)
SELECT gene, conf, category, grp2, count(*) AS n_genomes, sum(present::INT) AS n_present, avg(present::INT) AS frac
FROM u GROUP BY ALL;

CREATE OR REPLACE TABLE bfd_scratch.kx2_ku_ctrl_cov AS
SELECT m.locustag, m.genus, m.grp, m.set_name, sum(m.present::INT) AS ctrl_found, count(*) AS ctrl_n,
       any_value(m.busco_complete) AS busco_complete, any_value(m.busco_dup) AS busco_dup, any_value(m.n50) AS n50,
       any_value(m.n_proteins) AS n_proteins, any_value(c.frac_prot_pfam) AS frac_prot_pfam, any_value(c.panel_cov) AS panel_cov
FROM bfd_scratch.kx2_ku_matrix m JOIN bfd_scratch.kx2_ku_cov c USING (locustag)
WHERE m.is_ctrl GROUP BY m.locustag, m.genus, m.grp, m.set_name;

CREATE OR REPLACE TABLE bfd_scratch.kx2_ku_ctrl_grp AS
SELECT CASE WHEN genus='Coemansia' THEN 'Coemansia_strains' ELSE grp END AS grp2,
       count(*) AS n_genomes, avg(ctrl_found*1.0/ctrl_n) AS mean_ctrl_detect, 1 - avg(ctrl_found*1.0/ctrl_n) AS false_absence_est,
       avg(CASE WHEN ctrl_found < ctrl_n THEN 1.0 ELSE 0 END) AS frac_genomes_any_ctrl_missing
FROM bfd_scratch.kx2_ku_ctrl_cov
WHERE set_name='backbone' OR genus='Coemansia'
GROUP BY 1;
