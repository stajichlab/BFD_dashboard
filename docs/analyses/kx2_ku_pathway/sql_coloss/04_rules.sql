-- 04_rules.sql
-- Purpose : the gene-identification rule table (one row per gene). The rule text IS the documented method.
-- Inputs  : none (hand-coded). Every Pfam accession was read from the 'pfam' field of a named SwissProt reference in
--           bfd.swissprot_subject (e.g. Ku: YKU70 P32807, ku70 Q5AVC7, XRCC6 P12956; PRKDC P78527) and its presence
--           in bfd.pfam was checked in 02.
-- Output  : bfd_scratch.kx2_ku_rules
-- Columns : gene, category, mode (both = Pfam rule AND name rule; either = Pfam OR name; pfam; name),
--           req_all (all accessions needed), req_any (>=1 needed if non-empty), excl (none may be present),
--           ie_thr (domain i-evalue cutoff for the Pfam rule), name_re / name_excl_re (RE2 on
--           'gene_name | protein_name' of the SwissProt best hit), conf (high/medium/low), is_ctrl, note.
-- Thresholds: Pfam ie_thr 1e-5 (1e-3 for short/diverged families); name rules need best hit pident>=25, scov>=40,
--           -log10(e)>=10 (applied in 05).
-- Reports : the gene-list table in README.
CREATE OR REPLACE TABLE bfd_scratch.kx2_ku_rules AS
SELECT * FROM (VALUES
-- (a) classical NHEJ
 ('Ku','NHEJ','pfam',['PF02735','PF03731'],[]::VARCHAR[],[]::VARCHAR[],1e-5,NULL,NULL,'high',false,'Ku core + Ku_N; Ku70 and Ku80 not separated'),
 ('Ku70','NHEJ','both',['PF02735','PF03731'],[]::VARCHAR[],[]::VARCHAR[],1e-5,'(?i)ku70|yku70|mus-51|xrcc6|hdf1|helicase II subunit 1',NULL,'high',false,'Ku core + best hit named Ku70-type'),
 ('Ku80','NHEJ','both',['PF02735','PF03731'],[]::VARCHAR[],[]::VARCHAR[],1e-5,'(?i)ku80|yku80|mus-52|xrcc5|hdf2|helicase II subunit 2',NULL,'high',false,'Ku core + best hit named Ku80-type'),
 ('DNA-PKcs','NHEJ','pfam',['PF00454'],['PF08163','PF20500','PF20502','PF19704'],[]::VARCHAR[],1e-5,NULL,NULL,'high',false,'PIKK kinase + >=1 DNA-PKcs-specific domain (read from human PRKDC P78527)'),
 ('Lig4','NHEJ','both',['PF04675','PF01068'],[]::VARCHAR[],[]::VARCHAR[],1e-5,'(?i)DNA ligase 4|DNA ligase IV|lig4|dnl4|mus-53',NULL,'high',false,'ligase N+M and best hit ligase 4 (separates from Lig1/Lig3)'),
 ('XRCC4/Lif1','NHEJ','either',['PF06632'],[]::VARCHAR[],[]::VARCHAR[],1e-3,'(?i)xrcc4|lif1|non-homologous end-joining factor LIF1','(?i)like|xlf','low',false,'short coiled coil; fungal Lif1 has no Pfam model, name rule fills in'),
 ('XLF/Nej1','NHEJ','either',['PF09302'],[]::VARCHAR[],[]::VARCHAR[],1e-3,'(?i)xlf|xrcc4-like|nej1|cernunnos|NHEJ1|non-homologous end-joining protein 1',NULL,'low',false,'short, diverged'),
 ('PAXX','NHEJ','pfam',['PF15384'],[]::VARCHAR[],[]::VARCHAR[],1e-3,NULL,NULL,'low',false,'short, probably animal-restricted; expected absent almost everywhere'),
 ('SNM1/Pso2/Artemis family','NHEJ','pfam',['PF07522','PF00753'],[]::VARCHAR[],[]::VARCHAR[],1e-5,NULL,NULL,'medium',false,'DRMBL + metallo-beta-lactamase; family level, Artemis not separable by Pfam'),
 ('Artemis (DCLRE1C)','NHEJ','both',['PF07522'],[]::VARCHAR[],[]::VARCHAR[],1e-5,'(?i)artemis|DCLRE1C',NULL,'low',false,'DRMBL + best hit Artemis (vertebrate reference only)'),
 ('Pol4 (X-family Pol mu/lambda)','NHEJ','pfam',['PF14792','PF14791'],[]::VARCHAR[],[]::VARCHAR[],1e-5,NULL,NULL,'medium',false,'X-family palm+thumb; Pol beta/TdT-like also match'),
 ('PNKP/Tpp1','NHEJ','pfam',['PF08645'],[]::VARCHAR[],[]::VARCHAR[],1e-5,NULL,NULL,'medium',false,'polynucleotide kinase 3-prime phosphatase domain'),
 ('TDP1','NHEJ','pfam',['PF06087'],[]::VARCHAR[],[]::VARCHAR[],1e-5,NULL,NULL,'medium',false,'Tyr-DNA phosphodiesterase 1'),
 ('TDP2','NHEJ','both',['PF03372'],[]::VARCHAR[],[]::VARCHAR[],1e-5,'(?i)tyrosyl-DNA phosphodiesterase 2|tdp2|ttrap',NULL,'low',false,'generic Exo_endo_phos; needs name'),
 ('Aprataxin (Hnt3)','NHEJ','pfam',['PF01230','PF16278'],[]::VARCHAR[],[]::VARCHAR[],1e-5,NULL,NULL,'medium',false,'HIT + zf-C2HE'),
-- (b) resection / MRN
 ('Mre11','MRN','pfam',['PF00149','PF04152'],[]::VARCHAR[],[]::VARCHAR[],1e-5,NULL,NULL,'high',false,'nuclease + DNA-binding'),
 ('Rad50','MRN','either',['PF04423','PF13476'],[]::VARCHAR[],[]::VARCHAR[],1e-3,'(?i)^rad50|rad50 ',NULL,'high',false,'Zn hook + AAA_23, or best hit Rad50'),
 ('Nbs1/Xrs2','MRN','name',[]::VARCHAR[],[]::VARCHAR[],[]::VARCHAR[],1e-5,'(?i)nibrin|nbs1|xrs2',NULL,'low',false,'no distinctive Pfam in fungi (FHA/BRCT generic); name only'),
 ('Sae2/CtIP','MRN','pfam',['PF08573'],[]::VARCHAR[],[]::VARCHAR[],1e-3,NULL,NULL,'medium',false,'SAE2 domain'),
 ('Exo1','MRN','both',['PF00867','PF00752'],[]::VARCHAR[],[]::VARCHAR[],1e-5,'(?i)exonuclease 1|exodeoxyribonuclease 1|exo1','(?i)flap|fen1|XPG|GEN1|rad2|rad27','medium',false,'XPG_N+XPG_I and best hit Exo1 (separates from FEN1/XPG/GEN1)'),
 ('Dna2','MRN','pfam',['PF08696','PF13086'],[]::VARCHAR[],[]::VARCHAR[],1e-5,NULL,NULL,'high',false,'Dna2 nuclease + AAA_11'),
-- (c) apical kinases and checkpoint
 ('ATM/Tel1','KINASE','pfam',['PF00454','PF11640'],[]::VARCHAR[],[]::VARCHAR[],1e-5,NULL,NULL,'high',false,'PIKK + TAN'),
 ('ATR/Mec1','KINASE','pfam',['PF00454','PF08064'],[]::VARCHAR[],[]::VARCHAR[],1e-5,NULL,NULL,'high',false,'PIKK + UME'),
 ('ATRIP/Ddc2','KINASE','either',['PF09798'],[]::VARCHAR[],[]::VARCHAR[],1e-3,'(?i)atrip|ddc2|lcd1','(?i)CSB|helicase','low',false,'diverged; HEAT_LCD1 (PF09798) or name'),
 ('53BP1/Rad9/Crb2','KINASE','either',['PF08605'],[]::VARCHAR[],[]::VARCHAR[],1e-3,'(?i)DNA repair protein RAD9$|crb2|53-binding|TP53-binding',NULL,'low',false,'Rad9_Rad53_bind or name; 9-1-1 clamp RAD9A excluded by anchoring'),
 ('Rad53/Chk2','KINASE','both',['PF00498','PF00069'],[]::VARCHAR[],[]::VARCHAR[],1e-5,'(?i)rad53|chk2|cds1|checkpoint kinase 2',NULL,'medium',false,'FHA+kinase and best hit Rad53/Chk2/Cds1'),
-- (d) telomere
 ('TERT','TELO','pfam',['PF12009'],[]::VARCHAR[],[]::VARCHAR[],1e-5,NULL,NULL,'high',false,'Telomerase_RBD'),
 ('Est1','TELO','pfam',['PF10374','PF10373'],[]::VARCHAR[],['PF13638'],1e-5,NULL,NULL,'medium',false,'EST1 + EST1_DNA_bind, no PIN domain (SMG6 has PIN)'),
 ('SMG6 (comparator, NMD)','TELO','pfam',['PF10374','PF10373','PF13638'],[]::VARCHAR[],[]::VARCHAR[],1e-5,NULL,NULL,'medium',false,'EST1 + PIN; within-family comparator'),
 ('TPP1/Tpz1/Est3 (PF10341)','TELO','pfam',['PF10341'],[]::VARCHAR[],[]::VARCHAR[],1e-3,NULL,NULL,'low',false,'one OB-fold family covers TPP1, Tpz1, Est3'),
 ('POT1/Cdc13-OB (PF02765)','TELO','pfam',['PF02765'],[]::VARCHAR[],[]::VARCHAR[],1e-3,NULL,NULL,'low',false,'OB-fold; Pot1 and Cdc13 share it'),
 ('Pot1','TELO','both',['PF02765'],[]::VARCHAR[],[]::VARCHAR[],1e-3,'(?i)protection of telomeres protein 1|pot1','(?i)potassium|thiolase','low',false,'PF02765 + best hit Pot1'),
 ('Rap1','TELO','either',['PF08914'],['PF09197','PF11626','PF21411'],[]::VARCHAR[],1e-3,'(?i)DNA-binding protein rap1','(?i)gtpase|gap|gds','low',false,'Myb_Rap1 (+ others) or name'),
 ('Rif1','TELO','pfam',['PF12231'],[]::VARCHAR[],[]::VARCHAR[],1e-5,NULL,NULL,'medium',false,'Rif1_N'),
 ('Taz1/TRF','TELO','pfam',['PF08558'],[]::VARCHAR[],[]::VARCHAR[],1e-3,NULL,NULL,'low',false,'TRF homology domain'),
 ('Stn1','TELO','pfam',[]::VARCHAR[],['PF10451','PF09170'],[]::VARCHAR[],1e-3,NULL,NULL,'low',false,'CST'),
 ('Ten1','TELO','pfam',[]::VARCHAR[],['PF15490','PF12658'],[]::VARCHAR[],1e-3,NULL,NULL,'low',false,'CST; short'),
 ('Cdc13','TELO','name',[]::VARCHAR[],[]::VARCHAR[],[]::VARCHAR[],1e-5,'(?i)cell division control protein 13$','(?i)cyclin','low',false,'Cdc13 Pfam models have no hits in this panel; name only'),
 ('CTC1','TELO','pfam',[]::VARCHAR[],['PF15489','PF15491'],[]::VARCHAR[],1e-3,NULL,NULL,'low',false,'CST large subunit'),
-- ubiquitous controls
 ('CTRL Rad51','CTRL','both',['PF08423','PF14520'],[]::VARCHAR[],[]::VARCHAR[],1e-5,'(?i)rad51','(?i)rad51[a-d]|rad51ap|dmc1','high',true,'HR recombinase (Dmc1/Rad57 paralogs excluded by name)'),
 ('CTRL Pol delta (Pol3)','CTRL','both',['PF00136','PF03104','PF14260'],[]::VARCHAR[],[]::VARCHAR[],1e-5,'(?i)polymerase delta|pol3|POLD1',NULL,'high',true,'replicative polymerase delta catalytic subunit'),
 ('CTRL Pol epsilon (Pol2)','CTRL','pfam',['PF08490','PF22634'],[]::VARCHAR[],[]::VARCHAR[],1e-5,NULL,NULL,'high',true,'polymerase epsilon catalytic subunit'),
 ('CTRL Pol alpha (Pol1)','CTRL','pfam',['PF12254','PF00136'],[]::VARCHAR[],[]::VARCHAR[],1e-5,NULL,NULL,'high',true,'polymerase alpha catalytic subunit'),
 ('CTRL Mcm2','CTRL','pfam',['PF12619','PF00493'],[]::VARCHAR[],[]::VARCHAR[],1e-5,NULL,NULL,'high',true,'replicative helicase subunit'),
 ('CTRL Rad54','CTRL','pfam',['PF08658','PF00176'],[]::VARCHAR[],[]::VARCHAR[],1e-5,NULL,NULL,'high',true,'HR translocase'),
 ('CTRL Rad52','CTRL','pfam',['PF04098'],[]::VARCHAR[],[]::VARCHAR[],1e-5,NULL,NULL,'medium',true,'HR mediator (can be diverged)')
) t(gene, category, mode, req_all, req_any, excl, ie_thr, name_re, name_excl_re, conf, is_ctrl, note);
