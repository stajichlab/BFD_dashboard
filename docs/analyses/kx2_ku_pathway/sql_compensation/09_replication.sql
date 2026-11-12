-- Purpose: replication test. Every fungal genus without Ku (genus: <50% of q_loose genomes carry a Ku domain) is assigned to an independent loss lineage and
--          compared with Ku-retaining relatives of the same order (or phylum when < 3 reference genera in the order), feature by feature.
-- Inputs: kx2_kc_genus_long, kx2_kc_wide, kx2_kc_grp
-- Outputs: bfd_scratch.kx2_kc_kuneg (Ku-negative genera and lineage labels), bfd_scratch.kx2_kc_repl (lineage x feature: lineage value, reference median, difference, n reference genera),
--          bfd_scratch.kx2_kc_repl_sign (per feature: lineages higher / lower / equal vs relatives)
-- Lineage rule (manual, no tree for the full set): Kickxellales + Harpellales = ONE lineage (loss on shared stem in backbone tree; reference = Dimargaritales+Ramicandelaberales);
--          Melampsora, Amoeboaphelidium, Basidiobolus, Neoconidiobolus, Piptocephalis = one lineage each (genus-level; independence is assumed, not tested).
-- Thresholds: q_loose genomes; reference = Ku-positive genera (frac Ku>=0.5) in the same ORD if >= 3 such genera, else same PHYLUM.
--   (Kx+Hp has no Ku+ genus in its orders, so its reference is the Ku+ Kickxellomycota = Dimargaris, Dispira, Ramicandelaber.)
-- Produces: n independent Ku-negative lineages and the sign test per feature. Effective n is the number of lineages, not genomes.
CREATE OR REPLACE TABLE kx2_kc_kuneg AS
SELECT l.PHYLUM, l.GENUS, l.n_genomes, l.frac_pos AS frac_ku_pos, m.ORD,
 CASE WHEN m.ORD IN ('Kickxellales','Harpellales') THEN 'Kx+Hp' ELSE l.GENUS END AS lineage
FROM kx2_kc_genus_long l JOIN (SELECT PHYLUM, GENUS, any_value(ORD) AS ORD FROM kx2_kc_wide GROUP BY 1,2) m USING (PHYLUM, GENUS)
WHERE l.feature='ku' AND l.frac_pos<0.5;

CREATE OR REPLACE TABLE kx2_kc_repl AS
WITH feats AS (SELECT * FROM (VALUES ('polq'),('parp1l'),('parp_any'),('lig3'),('lig4'),('ligfam'),('rad51'),('rad51fam'),('rad51para'),('dmc1'),('rad52fam'),('rad54N'),('brca2dom'),('mus81'),('mms4'),('xpf'),('ercc4fam'),('fen1'),('tert'),('tertdom'),('n_genes')) t(feature)),
 gord AS (SELECT PHYLUM, GENUS, any_value(ORD) AS ORD FROM kx2_kc_wide GROUP BY 1,2),
 kupos AS (SELECT l.PHYLUM, l.GENUS, o.ORD FROM kx2_kc_genus_long l JOIN gord o USING (PHYLUM, GENUS) WHERE l.feature='ku' AND l.frac_pos>=0.5),
 lin AS (SELECT n.lineage, any_value(n.PHYLUM) AS PHYLUM, any_value(n.ORD) AS ORD FROM kx2_kc_kuneg n GROUP BY 1),
 scope AS (SELECT lin.lineage, lin.PHYLUM, lin.ORD, (SELECT count(*) FROM kupos k WHERE k.PHYLUM=lin.PHYLUM AND k.ORD=lin.ORD)>=3 AS use_ord FROM lin),
 refset AS (SELECT s.lineage, k.PHYLUM, k.GENUS FROM scope s JOIN kupos k ON k.PHYLUM=s.PHYLUM AND (NOT s.use_ord OR k.ORD=s.ORD)),
 lv AS (SELECT n.lineage, f.feature, median(g.med) AS lineage_val FROM kx2_kc_kuneg n JOIN feats f ON true JOIN kx2_kc_genus_long g ON g.PHYLUM=n.PHYLUM AND g.GENUS=n.GENUS AND g.feature=f.feature GROUP BY 1,2),
 rv AS (SELECT r.lineage, f.feature, median(g.med) AS ref_val, count(distinct g.GENUS) AS n_ref FROM refset r JOIN feats f ON true JOIN kx2_kc_genus_long g ON g.PHYLUM=r.PHYLUM AND g.GENUS=r.GENUS AND g.feature=f.feature GROUP BY 1,2)
SELECT lv.lineage, lv.feature, lv.lineage_val, rv.ref_val, lv.lineage_val-rv.ref_val AS diff, rv.n_ref FROM lv JOIN rv USING (lineage, feature);

CREATE OR REPLACE TABLE kx2_kc_repl_sign AS
SELECT feature, count(*) AS lineages, sum((diff>0)::int) AS higher, sum((diff<0)::int) AS lower, sum((diff=0)::int) AS equal FROM kx2_kc_repl GROUP BY 1;
SELECT * FROM kx2_kc_repl_sign ORDER BY 1;
