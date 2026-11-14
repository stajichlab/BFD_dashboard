-- 06_matrix.sql
-- Purpose : presence/absence matrix gene x genome (all 46 rules x 169 genomes) with QC covariates attached.
-- Inputs  : bfd_scratch.kx2_ku_rules, kx2_ku_genomes, kx2_ku_calls
-- Output  : bfd_scratch.kx2_ku_matrix (gene, category, conf, is_ctrl, locustag, genus, grp, set_name, present (bool),
--           n_present (proteins meeting the rule), n_pfam_only, n_name_only, QC covariates)
-- Thresholds: present = at least one protein with present=true in 05.
-- Reports : presence by group (07), all downstream tables.
CREATE OR REPLACE TABLE bfd_scratch.kx2_ku_matrix AS
SELECT r.gene, r.category, r.conf, r.is_ctrl, g.locustag, g.genus, g.ord, g.grp, g.set_name, g.role,
       coalesce(c.n_present,0) > 0 AS present,
       coalesce(c.n_present,0) AS n_present, coalesce(c.n_pfam_hit,0) AS n_pfam_hit, coalesce(c.n_name_hit,0) AS n_name_hit,
       g.busco_complete, g.busco_dup, g.n50, g.n_proteins
FROM bfd_scratch.kx2_ku_rules r CROSS JOIN bfd_scratch.kx2_ku_genomes g
LEFT JOIN (
  SELECT gene, locustag, count(*) FILTER (WHERE present) AS n_present,
         count(*) FILTER (WHERE pfam_ok) AS n_pfam_hit, count(*) FILTER (WHERE name_ok) AS n_name_hit
  FROM bfd_scratch.kx2_ku_calls GROUP BY 1,2
) c ON c.gene=r.gene AND c.locustag=g.locustag;
