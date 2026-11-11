-- Purpose: fusion probe (FAS, AROM, CAD) by Pfam co-occurrence within one protein.
-- Inputs: kx2_h1_prot. Output: kx2_h1_fusion. Report: fusion table (fraction of genomes per genus with fused form, median over genera).
-- Markers: FAS fused = PF16073(SAT)+PF00109(KS)+(PF18314|PF18325); alpha-only = KS+(PF18314|PF18325), no SAT; beta-only = SAT, no KS.
--          AROM modules PF01761, PF00275, PF01202, PF01487, PF01488|PF08501; fused = >=4 modules in one protein.
--          CAD fused = PF02786 + PF00185 in one protein; cps_only / atc_only otherwise.
-- Markers were chosen from the modal architecture of SwissProt-anchored FAS/AROM/URA2 hits (see 14_fusion_discovery.sql).
-- Filter: q genomes plus all Harpellales/Dimargaritales.
CREATE OR REPLACE TABLE kx2_h1_fusion AS
WITH D AS (SELECT lt, protein_id, length, list_distinct(string_split(arch,';')) AS ds FROM kx2_h1_prot WHERE arch IS NOT NULL),
F AS (SELECT lt, protein_id,
 list_contains(ds,'PF16073') AS sat, list_contains(ds,'PF00109') AS ks, (list_contains(ds,'PF18314') OR list_contains(ds,'PF18325')) AS alp,
 (list_contains(ds,'PF01761')::int+list_contains(ds,'PF00275')::int+list_contains(ds,'PF01202')::int+list_contains(ds,'PF01487')::int+(list_contains(ds,'PF01488') OR list_contains(ds,'PF08501'))::int) AS arom_n,
 list_contains(ds,'PF02786') AS cpl, list_contains(ds,'PF00185') AS otc, list_contains(ds,'PF00117') AS gat FROM D)
SELECT lt AS LOCUSTAG,
 count(*) FILTER (WHERE sat AND ks AND alp) AS fas_fused, count(*) FILTER (WHERE alp AND ks AND NOT sat) AS fas_alpha_only, count(*) FILTER (WHERE sat AND NOT ks) AS fas_beta_only,
 count(*) FILTER (WHERE arom_n>=4) AS arom_fused, count(*) FILTER (WHERE arom_n BETWEEN 1 AND 2) AS arom_module_split, count(*) FILTER (WHERE arom_n=3) AS arom_3,
 count(*) FILTER (WHERE cpl AND otc) AS cad_fused, count(*) FILTER (WHERE cpl AND NOT otc) AS cps_only, count(*) FILTER (WHERE otc AND NOT cpl) AS atc_only
FROM F GROUP BY 1;
WITH x AS (SELECT g.grp, g.GENUS, g.q, f.* FROM kx2_h1_fusion f JOIN kx2_h1_genomes g USING(LOCUSTAG) WHERE g.q OR g.grp IN ('Harpellales','Dimargaritales')),
y AS (SELECT grp, GENUS, count(*) n, avg((fas_fused>0)::int) fas_fused, avg((fas_alpha_only>0 AND fas_fused=0)::int) fas_split, avg((arom_fused>0)::int) arom_fused, avg((arom_fused=0 AND arom_module_split>0)::int) arom_split_only, avg((cad_fused>0)::int) cad_fused, avg((cad_fused=0 AND cps_only>0 AND atc_only>0)::int) cad_split FROM x GROUP BY 1,2)
SELECT grp, count(*) ngen, sum(n) ngenomes, round(median(fas_fused),2) fas_fused, round(median(fas_split),2) fas_split, round(median(arom_fused),2) arom_fused, round(median(arom_split_only),2) arom_split, round(median(cad_fused),2) cad_fused, round(median(cad_split),2) cad_split FROM y GROUP BY 1 ORDER BY 1;
