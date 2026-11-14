-- Purpose: artefact / robustness checks behind the report: (a) quality-strict sensitivity for Kickxellomycota, (b) Rad54 architecture (Rad54_N detection) in KxHp vs relatives,
--          (c) class composition of Rad51-domain proteins (Dmc1/Rad57 loss in Kickxellales), (d) Rad54-like presence vs BUSCO within Coemansia.
-- Inputs: kx2_kc_wide, kx2_kc_calls, kx2_kc_genomes
-- Outputs: bfd_scratch.kx2_kc_sens_strict, kx2_kc_rad54_arch, kx2_kc_rad51dom_cls, kx2_kc_rad54N_busco
-- Thresholds: q_loose (BUSCO >= 60, >= 3000 genes, dup <= 10) vs q_strict (BUSCO >= 85, dup <= 10); Rad54 class = best SwissProt hit 'DNA repair protein rhp54' (+ rad54/rdh54 names) AND SNF2-rel_dom or Rad54_N.
-- Produces: strict-subset prevalences (Kickxellales q_strict n=40: Ku 0, Lig3 0, Lig4 0, POLQ 0, Dmc1 0.03), Rad54_N detection 29 of 148 KxHp Rad54-like proteins vs 283/312 in Mucoromycota.
CREATE OR REPLACE TABLE kx2_kc_sens_strict AS
SELECT ORD, q_strict, count(*) AS n, avg(complete_pct) AS busco, avg((ku>0)::int) AS ku, avg((lig3>0)::int) AS lig3, avg((lig4>0)::int) AS lig4, avg((parp1l>0)::int) AS parp1l, avg((polq>0)::int) AS polq,
 avg(rad51fam) AS rad51fam, avg(rad51para) AS rad51para, avg((dmc1>0)::int) AS dmc1, avg((rad54N>0)::int) AS rad54N, avg((rad54name>0)::int) AS rad54name, avg(rad52fam) AS rad52fam, avg(tert) AS tert
FROM kx2_kc_wide WHERE ORD IN ('Kickxellales','Harpellales','Dimargaritales') AND q_loose GROUP BY ALL;

CREATE OR REPLACE TABLE kx2_kc_rad54_arch AS
SELECT CASE WHEN k.ORD IN ('Kickxellales','Harpellales') THEN 'KxHp' WHEN k.ORD IN ('Dimargaritales','Ramicandelaberales') THEN 'DimRam' WHEN k.PHYLUM='Mucoromycota' THEN 'Mucoromycota' ELSE 'other' END AS g,
 regexp_replace(coalesce(c.arch,'none'),'>?DEAD','') AS arch, count(*) AS n, count(distinct c.LOCUSTAG) AS genomes, avg(c.sp_lev) AS mean_lev, avg(c.qcov_x10)/10 AS mean_qcov
FROM kx2_kc_calls c JOIN kx2_kc_genomes k USING (LOCUSTAG)
WHERE c.cls IN ('RAD54','RDH54','RAD54B') AND c.sp_name ILIKE '%rhp54%' AND k.PHYLUM IN ('Kickxellomycota','Mucoromycota') GROUP BY ALL;

CREATE OR REPLACE TABLE kx2_kc_rad51dom_cls AS
WITH r AS (SELECT c.protein_id, c.LOCUSTAG, coalesce(c.cls,'(unclassified)') AS cls FROM kx2_kc_calls c WHERE list_contains(c.doms,'Rad51')),
 gg AS (SELECT LOCUSTAG, CASE WHEN ORD IN ('Kickxellales','Harpellales','Dimargaritales') THEN ORD WHEN PHYLUM IN ('Mucoromycota','Zoopagomycota') THEN PHYLUM END AS g FROM kx2_kc_wide WHERE q_loose)
SELECT gg.g, r.cls, count(*) AS n, count(distinct r.LOCUSTAG) AS genomes, (SELECT count(*) FROM gg g2 WHERE g2.g=gg.g) AS genomes_in_group
FROM gg JOIN r USING (LOCUSTAG) WHERE gg.g IS NOT NULL GROUP BY ALL;

CREATE OR REPLACE TABLE kx2_kc_rad54N_busco AS
SELECT (rad54N>0) AS has_rad54N, count(*) AS n, avg(complete_pct) AS busco, avg(n_genes) AS genes, avg(N50_bp) AS n50
FROM kx2_kc_wide WHERE GENUS='Coemansia' AND q_loose GROUP BY 1;
SELECT * FROM kx2_kc_rad54N_busco;
