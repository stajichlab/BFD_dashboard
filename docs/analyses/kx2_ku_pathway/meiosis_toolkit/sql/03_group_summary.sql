-- 03_group_summary.sql
-- Purpose : (a) control miss rates per group (single-domain and multi-domain style); (b) per gene and group: genomes with the gene / genomes
--           (annotation level); (c) an independent recount of four genes straight from bfd.pfam (verification); (d) Spo11 domain split.
-- Inputs  : bfd_scratch.kx2_mt_matrix (02), bfd_scratch.kx2_mt_pfam (extraction of bfd.pfam), bfd_scratch.kx2_ku_genomes, bfd.pfam.
-- Thresholds: a gene's absence in a group is INFORMATIVE only if the group's single-domain control miss rate is <= 10% (fixed in 01_rules.sql).
-- Reports : README tables 'controls', 'annotation presence by group', 'Spo11 domains'.
-- (a) controls
SELECT grp, count(DISTINCT locustag) AS n_genomes,
       round(100.0 * avg(CASE WHEN style='single' THEN CASE WHEN n_present=0 THEN 1.0 ELSE 0.0 END END), 1) AS single_miss_pct,
       round(100.0 * avg(CASE WHEN style='multi'  THEN CASE WHEN n_present=0 THEN 1.0 ELSE 0.0 END END), 1) AS multi_miss_pct
FROM bfd_scratch.kx2_mt_matrix WHERE is_ctrl GROUP BY grp ORDER BY grp;

-- (b) presence by group (present / genomes)
SELECT gene, style, conf,
  count(*) FILTER (WHERE grp='Kickxellales_core' AND n_present>0) || '/' || count(*) FILTER (WHERE grp='Kickxellales_core') AS core_kx,
  count(*) FILTER (WHERE grp='Mycoemilia_Spiromyces' AND n_present>0) || '/' || count(*) FILTER (WHERE grp='Mycoemilia_Spiromyces') AS myco_spiro,
  count(*) FILTER (WHERE grp='Ramicandelaber' AND n_present>0) || '/' || count(*) FILTER (WHERE grp='Ramicandelaber') AS rami,
  count(*) FILTER (WHERE grp='Harpellales' AND n_present>0) || '/' || count(*) FILTER (WHERE grp='Harpellales') AS harp,
  count(*) FILTER (WHERE grp='Dimargaritales' AND n_present>0) || '/' || count(*) FILTER (WHERE grp='Dimargaritales') AS dimarg,
  count(*) FILTER (WHERE grp='Zoopagales' AND n_present>0) || '/' || count(*) FILTER (WHERE grp='Zoopagales') AS zoop,
  count(*) FILTER (WHERE grp='Entomophthorales_Basidiobolus' AND n_present>0) || '/' || count(*) FILTER (WHERE grp='Entomophthorales_Basidiobolus') AS ento,
  count(*) FILTER (WHERE grp='Mucoromycota_s.l.' AND n_present>0) || '/' || count(*) FILTER (WHERE grp='Mucoromycota_s.l.') AS mucor,
  count(*) FILTER (WHERE grp='Root' AND n_present>0) || '/' || count(*) FILTER (WHERE grp='Root') AS root
FROM bfd_scratch.kx2_mt_matrix WHERE NOT is_ctrl GROUP BY gene, style, conf, fn ORDER BY fn, gene;

-- (c) independent recount from bfd.pfam of Spo11 (both domains, strict), Msh5 (MutS5_N + MutS_V or MutS_III), Hop2, Mnd1 (matches (b))
WITH g AS (SELECT locustag, grp FROM bfd_scratch.kx2_ku_genomes),
p AS (SELECT split_part(protein_id,'_',1) AS locustag, protein_id, split_part(pfam_acc,'.',1) AS acc, domain_i_evalue AS ie
      FROM bfd.pfam WHERE split_part(protein_id,'_',1) IN (SELECT locustag FROM bfd_scratch.kx2_ku_genomes)
        AND split_part(pfam_acc,'.',1) IN ('PF04406','PF21180','PF27035','PF00488','PF05192','PF07106','PF03962')),
spo AS (SELECT DISTINCT locustag FROM p WHERE acc IN ('PF04406','PF21180') AND ie<=1e-5 GROUP BY locustag, protein_id HAVING count(DISTINCT acc)=2),
m5 AS (SELECT DISTINCT locustag FROM p WHERE acc='PF27035' AND ie<=1e-3 AND protein_id IN (SELECT protein_id FROM p WHERE acc IN ('PF00488','PF05192') AND ie<=1e-3)),
h2 AS (SELECT DISTINCT locustag FROM p WHERE acc='PF07106' AND ie<=1e-3), m1 AS (SELECT DISTINCT locustag FROM p WHERE acc='PF03962' AND ie<=1e-3)
SELECT g.grp, count(*) AS n, count(*) FILTER (WHERE g.locustag IN (SELECT locustag FROM spo)) AS spo11_recount,
       count(*) FILTER (WHERE g.locustag IN (SELECT locustag FROM m5)) AS msh5_recount,
       count(*) FILTER (WHERE g.locustag IN (SELECT locustag FROM h2)) AS hop2_recount,
       count(*) FILTER (WHERE g.locustag IN (SELECT locustag FROM m1)) AS mnd1_recount
FROM g GROUP BY g.grp ORDER BY g.grp;

-- (d) Spo11 in core Kickxellales: which domain is detected, and what the Toprim-only proteins are (best SwissProt hit name)
SELECT g.genus, count(DISTINCT g.locustag) AS n_genomes,
       count(DISTINCT g.locustag) FILTER (WHERE EXISTS (SELECT 1 FROM bfd_scratch.kx2_mt_pfam p WHERE p.locustag=g.locustag AND p.acc='PF21180')) AS toprim_any,
       count(DISTINCT g.locustag) FILTER (WHERE EXISTS (SELECT 1 FROM bfd_scratch.kx2_mt_pfam p WHERE p.locustag=g.locustag AND p.acc='PF04406')) AS tp6a_any
FROM bfd_scratch.kx2_ku_genomes g WHERE g.grp='Kickxellales_core' GROUP BY g.genus ORDER BY n_genomes DESC;
WITH t AS (SELECT DISTINCT p.protein_id, p.locustag FROM bfd_scratch.kx2_mt_pfam p JOIN bfd_scratch.kx2_ku_genomes g USING (locustag)
           WHERE p.acc='PF21180' AND g.grp='Kickxellales_core'
             AND NOT EXISTS (SELECT 1 FROM bfd_scratch.kx2_mt_pfam q WHERE q.protein_id=p.protein_id AND q.acc='PF04406'))
SELECT coalesce(s.protein_name, '(no SwissProt hit)') AS best_hit_name, count(*) AS n_proteins, count(DISTINCT t.locustag) AS n_genomes,
       round(median(s.pident),1) AS median_pident, round(median(s.scov),1) AS median_scov
FROM t LEFT JOIN bfd_scratch.kx2_ku_sp s ON s.protein_id=t.protein_id GROUP BY 1 ORDER BY 2 DESC LIMIT 8;
