-- Purpose: independent re-computation check. Recount key Pfam-defined gene/family numbers straight from bfd.pfam (bypassing kx2_kc_pfhits/arch/calls/gene_counts/featmap)
--          for 313 genomes (Kickxellales, Harpellales, Dimargaritales, Ramicandelaberales, Rhizopus, Umbelopsis) and compare with kx2_kc_wide; recompute telomere detection
--          for N50>=500 kb genomes in Kickxellomycota vs Mucoromycota+Mortierellomycota+Glomeromycota from bfd.telomere_summary and compare with kx2_kc_tel_genome.
-- Inputs: bfd.pfam, bfd.telomere_summary, kx2_kc_genomes, kx2_kc_wide, kx2_kc_tel_genome
-- Outputs: temporary bfd_scratch.kx2_kc_chk_pf, kx2_kc_chk_tel (dropped at the end)
-- Expected result (2026-11-11 run): 0 mismatches in 313 genomes for ku, rad51fam, rad54N, ligfam, parp_any, tertdom, rad52fam, polq; telomere detection identical
--          (Harpellales 3/3, Kickxellales 0/2, Mucoromycota+Mortierellomycota+Glomeromycota 12/63).
-- Not re-checked here: name-based classes (03/04) and the genus-level / Cliff tables (these were rebuilt once after the rhp54 name fix, not independently recomputed).
CREATE OR REPLACE TABLE kx2_kc_chk_pf AS
WITH g AS (SELECT LOCUSTAG FROM kx2_kc_genomes WHERE ORD IN ('Kickxellales','Harpellales','Dimargaritales','Ramicandelaberales') OR GENUS IN ('Rhizopus','Umbelopsis')),
p AS (SELECT split_part(protein_id,'_',1) AS LOCUSTAG, regexp_replace(protein_id,'-T[0-9]+$','') AS gene, pfam_id FROM bfd.pfam WHERE domain_i_evalue<=1e-5 AND pfam_id IN ('Ku','Rad51','Rad54_N','DNA_ligase_A_M','PARP','Telomerase_RBD','Rad52_Rad22','DNA_pol_A','DEAD','Helicase_C','POLQ_helical') AND split_part(protein_id,'_',1) IN (SELECT LOCUSTAG FROM g)),
pg AS (SELECT LOCUSTAG, gene, list(DISTINCT pfam_id) d FROM p GROUP BY 1,2)
SELECT g.LOCUSTAG,
 count(DISTINCT CASE WHEN list_contains(d,'Ku') THEN gene END) ku,
 count(DISTINCT CASE WHEN list_contains(d,'Rad51') THEN gene END) rad51fam,
 count(DISTINCT CASE WHEN list_contains(d,'Rad54_N') THEN gene END) rad54N,
 count(DISTINCT CASE WHEN list_contains(d,'DNA_ligase_A_M') THEN gene END) ligfam,
 count(DISTINCT CASE WHEN list_contains(d,'PARP') THEN gene END) parp_any,
 count(DISTINCT CASE WHEN list_contains(d,'Telomerase_RBD') THEN gene END) tertdom,
 count(DISTINCT CASE WHEN list_contains(d,'Rad52_Rad22') THEN gene END) rad52fam,
 count(DISTINCT CASE WHEN list_contains(d,'DNA_pol_A') AND (list_contains(d,'DEAD') OR list_contains(d,'Helicase_C') OR list_contains(d,'POLQ_helical')) THEN gene END) polq
FROM g LEFT JOIN pg USING (LOCUSTAG) GROUP BY 1;
CREATE OR REPLACE TABLE kx2_kc_chk_tel AS
SELECT CASE WHEN w.PHYLUM='Kickxellomycota' THEN w.ORD WHEN w.PHYLUM IN ('Mucoromycota','Mortierellomycota','Glomeromycota') THEN 'Mucoro' ELSE 'x' END grp, count(*) n, sum((t.telomere_tracts>0)::int) det
FROM kx2_kc_wide w LEFT JOIN bfd.telomere_summary t USING (LOCUSTAG) WHERE w.q_loose AND w.N50_bp>=500000 AND (w.PHYLUM IN ('Kickxellomycota','Mucoromycota','Mortierellomycota','Glomeromycota')) GROUP BY 1;
-- comparison (all mism must be 0)
SELECT 'ku' f, sum((c.ku<>w.ku)::int) mism, count(*) n FROM kx2_kc_chk_pf c JOIN kx2_kc_wide w USING (LOCUSTAG)
UNION ALL SELECT 'rad51fam', sum((c.rad51fam<>w.rad51fam)::int), count(*) FROM kx2_kc_chk_pf c JOIN kx2_kc_wide w USING (LOCUSTAG)
UNION ALL SELECT 'rad54N', sum((c.rad54N<>w.rad54N)::int), count(*) FROM kx2_kc_chk_pf c JOIN kx2_kc_wide w USING (LOCUSTAG)
UNION ALL SELECT 'ligfam', sum((c.ligfam<>w.ligfam)::int), count(*) FROM kx2_kc_chk_pf c JOIN kx2_kc_wide w USING (LOCUSTAG)
UNION ALL SELECT 'parp_any', sum((c.parp_any<>w.parp_any)::int), count(*) FROM kx2_kc_chk_pf c JOIN kx2_kc_wide w USING (LOCUSTAG)
UNION ALL SELECT 'tertdom', sum((c.tertdom<>w.tertdom)::int), count(*) FROM kx2_kc_chk_pf c JOIN kx2_kc_wide w USING (LOCUSTAG)
UNION ALL SELECT 'rad52fam', sum((c.rad52fam<>w.rad52fam)::int), count(*) FROM kx2_kc_chk_pf c JOIN kx2_kc_wide w USING (LOCUSTAG)
UNION ALL SELECT 'polq', sum((c.polq<>w.polq)::int), count(*) FROM kx2_kc_chk_pf c JOIN kx2_kc_wide w USING (LOCUSTAG);
SELECT c.grp, c.n chk_n, c.det chk_det, s.n stored_n, s.det stored_det FROM kx2_kc_chk_tel c LEFT JOIN (SELECT CASE WHEN PHYLUM='Kickxellomycota' THEN ORD ELSE 'Mucoro' END grp, count(*) n, sum(has_tel::int) det FROM kx2_kc_tel_genome WHERE mid_contig AND PHYLUM IN ('Kickxellomycota','Mucoromycota','Mortierellomycota','Glomeromycota') GROUP BY 1) s USING (grp) ORDER BY 1;
DROP TABLE kx2_kc_chk_pf;
DROP TABLE kx2_kc_chk_tel;
