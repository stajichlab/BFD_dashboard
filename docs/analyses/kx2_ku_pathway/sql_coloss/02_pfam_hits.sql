-- 02_pfam_hits.sql
-- Purpose : per-protein Pfam hits for the accession panel used by the Ku-pathway gene rules,
--           restricted to the 169 genomes of kx2_ku_genomes.
-- Inputs  : bfd.pfam (pfam_acc carries a version suffix, e.g. PF02735.19; stripped here),
--           bfd_scratch.kx2_ku_genomes
-- Output  : bfd_scratch.kx2_ku_accs  (panel of accessions; every accession was read from
--                                     bfd.swissprot_annot/subject 'pfam' field for a named reference protein)
--           bfd_scratch.kx2_ku_pfam   (protein_id, locustag, acc, pfam_id, min_ie = best domain i-evalue, n_dom,
--                                      first_from/last_to = envelope span on the protein)
-- Thresholds: keep rows with domain_i_evalue <= 1e-3 (gene rules then apply 1e-5 or 1e-3 per rule, see 03).
-- Reports : none directly (intermediate).
CREATE OR REPLACE TABLE bfd_scratch.kx2_ku_accs AS
SELECT * FROM (VALUES
 ('PF02735','Ku'),('PF03731','Ku_N'),('PF03730','Ku_C'),('PF08785','Ku_PK_bind'),('PF02037','SAP'),
 ('PF20500','DNA-PKcs_N'),('PF20502','DNAPKcs_CC1-2'),('PF08163','DNAPKcs_CC3'),('PF19704','DNAPKcs_CC5'),
 ('PF00454','PI3_PI4_kinase'),('PF02259','FAT'),('PF02260','FATC'),('PF11640','TAN'),('PF08064','UME'),
 ('PF23593','HEAT_ATR'),('PF25030','M-HEAT_ATR'),('PF25385','HEAT_MEC1_N'),
 ('PF04675','DNA_ligase_A_N'),('PF01068','DNA_ligase_A_M'),('PF04679','DNA_ligase_A_C'),('PF16589','BRCT_2'),('PF00533','BRCT'),('PF11411','DNA_ligase_IV'),
 ('PF06632','XRCC4'),('PF09302','XLF'),('PF15384','PAXX'),
 ('PF07522','DRMBL'),('PF00753','Lactamase_B'),
 ('PF14792','DNA_pol_B_palm'),('PF14791','DNA_pol_B_thumb'),('PF10391','DNA_pol_lambd_f'),('PF14716','HHH_8'),
 ('PF08645','PNK3P'),('PF13671','AAA_33'),('PF06087','Tyr-DNA_phospho'),('PF03372','Exo_endo_phos'),
 ('PF01230','HIT'),('PF16278','zf-C2HE'),('PF11969','DcpS_C'),
 ('PF00149','Metallophos'),('PF04152','Mre11_DNA_bind'),('PF13476','AAA_23'),('PF04423','Rad50_zn_hook'),('PF13558','SbcC_Walker_B'),
 ('PF00498','FHA'),('PF08599','Nbs1_C'),('PF16508','NIBRIN_BRCT_II'),('PF08573','SAE2'),
 ('PF00867','XPG_I'),('PF00752','XPG_N'),('PF08696','Dna2'),('PF13086','AAA_11'),('PF13087','AAA_12'),
 ('PF00069','Pkinase'),('PF08605','Rad9_Rad53_bind'),('PF22520','Crb2_Tudor'),('PF22395','Crb2-like_BRCT-like'),('PF09038','53-BP1_Tudor'),('PF18428','BRCT_3'),
 ('PF28426','ATRIP_C'),('PF28409','ATRIP_N'),('PF09798','LCD1'),
 ('PF12009','Telomerase_RBD'),('PF00078','RVT_1'),('PF11474','TEN_TERT'),('PF21399','TERT_C'),
 ('PF10374','EST1'),('PF10373','EST1_DNA_bind'),('PF13638','PIN_4'),('PF10341','TPP1'),
 ('PF02765','POT1'),('PF16686','POT1PC'),
 ('PF09197','Rap1-DNA-bind'),('PF11626','Rap1_C'),('PF08914','Myb_Rap1'),('PF21411','Rap1_C_Sp'),('PF00249','Myb_DNA-binding'),
 ('PF12231','Rif1_N'),('PF08558','TRF'),
 ('PF10451','Stn1'),('PF12659','Stn1_C'),('PF09170','STN1_2'),('PF12658','Ten1'),('PF15490','Ten1_2'),
 ('PF16853','CDC13_N'),('PF18691','Cdc13_OB2'),('PF18233','Cdc13_OB4_dimer'),('PF15489','CTC1'),('PF15491','CTC1_2'),
 -- controls
 ('PF08423','Rad51'),('PF14520','HHH_5'),
 ('PF00136','DNA_pol_B'),('PF03104','DNA_pol_B_exo1'),('PF14260','zf-C4pol'),('PF24055','POL3_N'),
 ('PF08490','DUF1744'),('PF22634','POL2_thumb'),('PF12254','DNA_pol_alpha_N'),
 ('PF00493','MCM'),('PF12619','MCM2_N'),('PF04098','Rad52_Rad22'),('PF08658','Rad54_N'),('PF00271','Helicase_C'),('PF00176','SNF2-rel_dom')
) t(acc, name);

CREATE OR REPLACE TABLE bfd_scratch.kx2_ku_pfam AS
SELECT p.protein_id, split_part(p.protein_id,'_',1) AS locustag,
       split_part(p.pfam_acc,'.',1) AS acc, any_value(p.pfam_id) AS pfam_id,
       min(p.domain_i_evalue) AS min_ie, count(*) AS n_dom,
       min(p.env_from) AS first_from, max(p.env_to) AS last_to
FROM bfd.pfam p
JOIN bfd_scratch.kx2_ku_accs a ON a.acc = split_part(p.pfam_acc,'.',1)
WHERE split_part(p.protein_id,'_',1) IN (SELECT locustag FROM bfd_scratch.kx2_ku_genomes)
  AND p.domain_i_evalue <= 1e-3
GROUP BY 1,2,3;
