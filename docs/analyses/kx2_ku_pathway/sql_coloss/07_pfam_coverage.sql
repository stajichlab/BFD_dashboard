-- 07_pfam_coverage.sql
-- Purpose : per-genome annotation-completeness covariates.
--           (a) frac_prot_pfam = share of proteins with any Pfam hit (domain i-evalue <= 1e-5);
--           (b) panel_cov     = share of a near-universal Pfam panel present in the genome. Panel = Pfam families found
--               (i-evalue <= 1e-5) in >= 95% of the 48 backbone genomes OUTSIDE Kickxellomycota
--               (Mucoromycota s.l., Zoopagales, Entomophthorales, Basidiobolus, root). Panel size is reported in README.
--               The panel is defined from non-ingroup genomes so Kickxellomycota gene loss cannot shape it, but
--               it is data-derived (not the curated ~1,400 list) - caveat stated in README.
-- Inputs  : bfd.pfam, bfd.gene_proteins (via kx2_ku_genomes.n_proteins), kx2_ku_genomes
-- Output  : bfd_scratch.kx2_ku_pfam_all (locustag, acc, n_prot), bfd_scratch.kx2_ku_panel (acc), bfd_scratch.kx2_ku_cov
-- Thresholds: 1e-5; panel prevalence >= 0.95 in 48 outgroup genomes.
-- Reports : QC covariates for the false-absence analysis (11).
CREATE OR REPLACE TABLE bfd_scratch.kx2_ku_pfam_all AS
SELECT split_part(protein_id,'_',1) AS locustag, split_part(pfam_acc,'.',1) AS acc, count(DISTINCT protein_id) AS n_prot
FROM bfd.pfam
WHERE split_part(protein_id,'_',1) IN (SELECT locustag FROM bfd_scratch.kx2_ku_genomes) AND domain_i_evalue <= 1e-5
GROUP BY 1,2;

CREATE OR REPLACE TABLE bfd_scratch.kx2_ku_panel AS
SELECT acc FROM bfd_scratch.kx2_ku_pfam_all a JOIN bfd_scratch.kx2_ku_genomes g USING (locustag)
WHERE g.set_name='backbone' AND g.phylum <> 'Kickxellomycota'
GROUP BY acc HAVING count(DISTINCT locustag) >= 0.95 * (SELECT count(*) FROM bfd_scratch.kx2_ku_genomes WHERE set_name='backbone' AND phylum <> 'Kickxellomycota');

CREATE OR REPLACE TABLE bfd_scratch.kx2_ku_cov AS
WITH pp AS (
  SELECT split_part(protein_id,'_',1) AS locustag, count(DISTINCT protein_id) AS n_prot_pfam
  FROM bfd.pfam WHERE split_part(protein_id,'_',1) IN (SELECT locustag FROM bfd_scratch.kx2_ku_genomes) AND domain_i_evalue <= 1e-5
  GROUP BY 1
),
pc AS (
  SELECT a.locustag, count(*) AS panel_hit FROM bfd_scratch.kx2_ku_pfam_all a JOIN bfd_scratch.kx2_ku_panel p USING (acc) GROUP BY 1
)
SELECT g.locustag, g.grp, g.set_name, g.genus, g.n_proteins, pp.n_prot_pfam, pp.n_prot_pfam * 1.0 / g.n_proteins AS frac_prot_pfam,
       pc.panel_hit, pc.panel_hit * 1.0 / (SELECT count(*) FROM bfd_scratch.kx2_ku_panel) AS panel_cov
FROM bfd_scratch.kx2_ku_genomes g LEFT JOIN pp USING (locustag) LEFT JOIN pc USING (locustag);
