-- Purpose: per-protein gene calls = best-hit name class (03) AND the required Pfam domain (02); plus Pfam-only architecture flags.
-- Inputs: kx2_kc_arch, kx2_kc_sbclass
-- Output: bfd_scratch.kx2_kc_calls (protein_id, arch, doms, ku_pf, polq_pf, parp_pf, parp1like_pf, cls, sp_name, pf_ok)
-- Rules: ku_pf = Ku domain or (Ku_N and Ku_C). polq_pf = DNA_pol_A AND (DEAD|Helicase_C|POLQ_helical) in the SAME protein (POLQ architecture).
--        parp1like_pf = PARP catalytic AND (WGR|zf-PARP|PARP_reg). pf_ok = class-specific required domain (see CASE).
-- Produces: POLQ architecture counts (27 proteins kingdom-wide), PARP1-like architecture, class calls.
-- Note: the name class 'POLQ' (best hit 'DNA polymerase theta') is mostly HELQ/Mus308-type helicases: only 5 of 6763 also carry DNA_pol_A, so POLQ is called from ARCHITECTURE (polq_pf), not from the name.
CREATE OR REPLACE TABLE kx2_kc_calls AS
WITH pf AS (
 SELECT protein_id, LOCUSTAG, arch, doms,
  (list_contains(doms,'Ku') OR (list_contains(doms,'Ku_N') AND list_contains(doms,'Ku_C'))) AS ku_pf,
  (list_contains(doms,'DNA_pol_A') AND (list_contains(doms,'DEAD') OR list_contains(doms,'Helicase_C') OR list_contains(doms,'POLQ_helical'))) AS polq_pf,
  list_contains(doms,'PARP') AS parp_pf,
  (list_contains(doms,'PARP') AND (list_contains(doms,'WGR') OR list_contains(doms,'zf-PARP') OR list_contains(doms,'PARP_reg'))) AS parp1like_pf
 FROM kx2_kc_arch)
SELECT coalesce(p.protein_id,s.protein_id) AS protein_id, coalesce(p.LOCUSTAG,s.LOCUSTAG) AS LOCUSTAG, p.arch, p.doms, p.ku_pf, p.polq_pf, p.parp_pf, p.parp1like_pf,
 s.cls, s.protein_name AS sp_name, s.neg_log10_evalue AS sp_lev, s.qcov_x10, s.scov_x10,
 CASE s.cls
  WHEN 'KU' THEN list_contains(p.doms,'Ku')
  WHEN 'POLQ' THEN list_contains(p.doms,'DNA_pol_A')
  WHEN 'PARP1' THEN list_contains(p.doms,'PARP')
  WHEN 'XRCC1' THEN list_contains(p.doms,'XRCC1_N')
  WHEN 'LIG1' THEN list_contains(p.doms,'DNA_ligase_A_M')
  WHEN 'LIG3' THEN list_contains(p.doms,'DNA_ligase_A_M')
  WHEN 'LIG4' THEN list_contains(p.doms,'DNA_ligase_A_M')
  WHEN 'FEN1' THEN list_contains(p.doms,'XPG_I')
  WHEN 'EXO1' THEN list_contains(p.doms,'XPG_I')
  WHEN 'GEN1' THEN list_contains(p.doms,'XPG_I')
  WHEN 'RAD51' THEN list_contains(p.doms,'Rad51')
  WHEN 'DMC1' THEN list_contains(p.doms,'Rad51')
  WHEN 'RAD51_PARALOG_BCD' THEN list_contains(p.doms,'Rad51')
  WHEN 'RAD55' THEN list_contains(p.doms,'Rad51') OR list_contains(p.doms,'ATPase')
  WHEN 'RAD57' THEN list_contains(p.doms,'Rad51') OR list_contains(p.doms,'ATPase')
  WHEN 'XRCC2' THEN list_contains(p.doms,'Rad51') OR list_contains(p.doms,'ATPase')
  WHEN 'XRCC3' THEN list_contains(p.doms,'Rad51') OR list_contains(p.doms,'ATPase')
  WHEN 'RAD52' THEN list_contains(p.doms,'Rad52_Rad22')
  WHEN 'RAD59' THEN list_contains(p.doms,'Rad52_Rad22')
  WHEN 'RAD54' THEN list_contains(p.doms,'SNF2-rel_dom') OR list_contains(p.doms,'Rad54_N')
  WHEN 'RAD54B' THEN list_contains(p.doms,'SNF2-rel_dom') OR list_contains(p.doms,'Rad54_N')
  WHEN 'RDH54' THEN list_contains(p.doms,'SNF2-rel_dom') OR list_contains(p.doms,'Rad54_N')
  WHEN 'BRCA2' THEN list_contains(p.doms,'BRCA2') OR list_contains(p.doms,'BRCA-2_OB1')
  WHEN 'MUS81' THEN list_contains(p.doms,'ERCC4')
  WHEN 'MMS4_EME1' THEN list_contains(p.doms,'ERCC4')
  WHEN 'XPF_RAD1' THEN list_contains(p.doms,'ERCC4')
  WHEN 'RAD10_ERCC1' THEN list_contains(p.doms,'Rad10')
  WHEN 'TERT' THEN list_contains(p.doms,'Telomerase_RBD') OR list_contains(p.doms,'RVT_1')
 END AS pf_ok
FROM pf p FULL JOIN kx2_kc_sbclass s ON s.protein_id=p.protein_id;
SELECT cls, sum(pf_ok::int) ok, count(*) FROM kx2_kc_calls WHERE cls IS NOT NULL GROUP BY 1 ORDER BY 1;
