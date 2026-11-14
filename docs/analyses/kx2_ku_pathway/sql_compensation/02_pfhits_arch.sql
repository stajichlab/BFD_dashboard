-- Purpose: pull Pfam domain hits for the DNA-repair anchor/architecture domains and build per-protein architecture strings.
-- Inputs: bfd.pfam (Pfam names verified against swissprot_subject.pfam of reference proteins, see README gene table)
-- Outputs: bfd_scratch.kx2_kc_pfhits (domain hits, domain_i_evalue<=1e-5), bfd_scratch.kx2_kc_arch (per protein: arch string ordered by ali_from, doms list)
-- Thresholds: domain_i_evalue <= 1e-5. HHH_5 kept in pfhits, dropped from arch strings.
-- Produces: architecture of POLQ, PARP, ligase, Rad54 etc. Pfam ids are matched on pfam_id (name); accession = split_part(pfam_acc,'.',1).
CREATE OR REPLACE TABLE kx2_kc_pfhits AS
WITH rel AS (SELECT * FROM (VALUES ('Ku'),('Ku_N'),('Ku_C'),('DEAD'),('Helicase_C'),('SNF2-rel_dom'),('Rad54_N'),('DNA_pol_A'),('POLQ_helical'),('PARP'),('WGR'),('zf-PARP'),('PARP_reg'),('BRCT'),('XRCC1_N'),('DNA_ligase_A_N'),('DNA_ligase_A_M'),('DNA_ligase_A_C'),('Rad51'),('Rad52_Rad22'),('ERCC4'),('Rad10'),('BRCA2'),('BRCA-2_OB1'),('XPG_N'),('XPG_I'),('Telomerase_RBD'),('RVT_1'),('TEN_TERT'),('HHH_5'),('ATPase')) t(pid))
SELECT p.protein_id, split_part(p.protein_id,'_',1) AS LOCUSTAG, p.pfam_id, split_part(p.pfam_acc,'.',1) AS acc, p.ali_from, p.ali_to, p.domain_i_evalue
FROM bfd.pfam p JOIN rel ON rel.pid=p.pfam_id WHERE p.domain_i_evalue<=1e-5;

CREATE OR REPLACE TABLE kx2_kc_arch AS
SELECT protein_id, LOCUSTAG, string_agg(pfam_id, '>' ORDER BY ali_from) AS arch, list(DISTINCT pfam_id) AS doms
FROM kx2_kc_pfhits WHERE pfam_id NOT IN ('HHH_5') GROUP BY 1,2;
-- check: accessions per domain
SELECT pfam_id, acc, count(*), count(distinct protein_id) FROM kx2_kc_pfhits GROUP BY ALL ORDER BY 3 DESC;
