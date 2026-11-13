-- 01_rules.sql  (written and fixed BEFORE any group-level result was viewed; 2026-11-12)
-- Purpose : gene-identification rules for the meiotic toolkit and its controls, one row per gene. The rule text IS the method.
-- Inputs  : none (hand-coded). Every Pfam accession and name below was read from the 'pfam' / 'protein_name' fields of named
--           SwissProt references in bfd.swissprot_subject (S. cerevisiae, S. pombe, human, mouse, A. thaliana, C. albicans), and
--           each accession's presence in bfd.pfam was confirmed (extraction in 02). References, for example:
--           Spo11 P23179 (TP6A_N PF04406 + TOP6A-Spo11_Toprim PF21180), Msh5 Q12175 (MutS5_N PF27035), Hop2 P53187 (WHD_TBPIP PF07106),
--           Mnd1 P53102 (Mnd1 PF03962), Mei4 P29467 (PF13971), Mei5 P32489 (PF10376), Mer2/Rec107 P21651 (PF09074), Rec114 P32841
--           (PF03525), Hop1 P20050 (HORMA PF02301), Red1 P14291 (PF07964), Pch2 P38126 (PF00004 + PF23242), Ndt80 P38830 (PF05224),
--           Spo16 P17122 (PF21698, PF19225), Mer3 P51979 (DEAD PF00270 + Helicase_C PF00271 + Sec63 PF02889).
-- Output  : bfd_scratch.kx2_mt_rules
-- Columns : gene, grp (function), mode (pfam | name | both | either), req_all (all accessions needed), req_any (>=1 needed if
--           non-empty), excl (none may be present), ie_thr (domain i-evalue cutoff; strict 1e-5 or loose 1e-3), name_re / name_excl_re
--           (RE2 on 'gene_name | protein_name' of the SwissProt best hit), conf, is_ctrl, style (single | multi domain control), note.
-- NOT ASSESSED (no usable rule): Zip1, Zip2, Zip3 (yeast proteins have no Pfam; the gene names collide with zinc transporters, ZIP
--           family), Msh4 / Mlh3 / Rec8 / Hop1 identity (close paralogs of mitotic genes, classified by trees in 05-06).
-- Thresholds fixed in advance (README): a meiotic gene's absence in a group is INFORMATIVE only if the matching single-domain control
--           (style 'single') is missed in <= 10% of that group's genomes.
-- REVISION 2026-11-12 (made after inspecting control miss rates and row counts only; NO meiotic-gene group result had been viewed):
--   * CTRL Mlh1: the first rule used MutL_C PF08676, which belongs to the Pms1/Mlh3 side; Mlh1 carries Mlh1_C PF16413 (read from the
--     human, yeast, S. pombe and Arabidopsis MLH1 references), so Pfam and name evidence never overlapped (0 present everywhere).
--   * CTRL Msh2 (name+MutS) dropped: only 8 best hits in the 169 genomes carry an Msh2 name, so the name evidence is too sparse; replaced
--     by a Pfam control 'CTRL MutS_I+II (Msh2/3/6)'.
--   * CTRL Rad50 (single) dropped: the Rad50 zinc hook model missed in 70-92% of genomes in EVERY group, a defective model, not a
--     measure of detection.
--   * Hop1 rule dropped: its name rule also matches the unrelated Arabidopsis protein named HOP1 (Hsp70-Hsp90 organizing protein), so
--     it cannot be called from annotation; Hop1 is handled by the HORMA family count and the tree.
CREATE OR REPLACE TABLE bfd_scratch.kx2_mt_rules AS
SELECT * FROM (VALUES
-- (A) double-strand-break formation
 ('Spo11','DSB','pfam',['PF04406','PF21180'],[]::VARCHAR[],[]::VARCHAR[],1e-5,NULL,NULL,'high',false,'multi','TopoVI-A-like N-terminus + Toprim; in fungi Spo11 is the only TopoVI-A-like protein (SPO11_like PF03533 is animal-specific)'),
 ('Rec114','DSB','pfam',[]::VARCHAR[],['PF03525','PF15165'],[]::VARCHAR[],1e-3,NULL,NULL,'low',false,'single','fast-evolving; yeast and animal Pfam models differ'),
 ('Mer2/Rec107','DSB','pfam',['PF09074'],[]::VARCHAR[],[]::VARCHAR[],1e-3,NULL,NULL,'low',false,'single','fast-evolving'),
 ('Mei4','DSB','pfam',['PF13971'],[]::VARCHAR[],[]::VARCHAR[],1e-3,NULL,NULL,'low',false,'single','S. pombe mei4 is an unrelated forkhead protein; Pfam model is yeast/animal Mei4'),
-- (B) strand exchange and cofactors (Dmc1 itself is in ../family_trees)
 ('Hop2','STRAND','pfam',['PF07106'],[]::VARCHAR[],[]::VARCHAR[],1e-3,NULL,NULL,'medium',false,'single','WHD_TBPIP'),
 ('Mnd1','STRAND','pfam',['PF03962'],[]::VARCHAR[],[]::VARCHAR[],1e-3,NULL,NULL,'medium',false,'single','Mnd1 (LZ3wCH PF18517 is shared with Hop2-like proteins and is not used)'),
 ('Mei5','STRAND','pfam',['PF10376'],[]::VARCHAR[],[]::VARCHAR[],1e-3,NULL,NULL,'low',false,'single','short, fast-evolving'),
 ('Swi5/Sae3','STRAND','pfam',['PF07061'],[]::VARCHAR[],[]::VARCHAR[],1e-3,NULL,NULL,'low',false,'single','Swi5 also acts in mitotic Rad51 mediation, so it is not meiosis-specific'),
-- (C) crossover / ZMM pathway
 ('Msh5','ZMM','pfam',['PF27035'],['PF00488','PF05192'],[]::VARCHAR[],1e-3,NULL,NULL,'medium',false,'multi','MutS5_N is found in Msh5 references only'),
 ('Mer3/Hfm1','ZMM','both',['PF00270','PF00271','PF02889'],[]::VARCHAR[],[]::VARCHAR[],1e-5,'(?i)mer3|hfm1','(?i)snrnp200|brr2|u5 small|cointegrator|polq|helicase q|ski2','medium',false,'multi','Sec63-type helicase; name rule separates it from Brr2/SNRNP200'),
 ('Spo16','ZMM','pfam',[]::VARCHAR[],['PF21698','PF19225','PF15162'],[]::VARCHAR[],1e-3,NULL,NULL,'low',false,'single','short, fast-evolving'),
-- (D) axis, checkpoint, cohesion
 ('Red1','AXIS','pfam',['PF07964'],[]::VARCHAR[],[]::VARCHAR[],1e-3,NULL,NULL,'low',false,'single','yeast/Schizosaccharomyces axis protein; fast-evolving'),
 ('Mek1','AXIS','both',['PF00498','PF00069'],[]::VARCHAR[],[]::VARCHAR[],1e-5,'(?i)mek1','(?i)rad53|chk2|cds1|checkpoint kinase','low',false,'multi','FHA + kinase, same architecture as Rad53'),
 ('Pch2/TRIP13','AXIS','pfam',['PF00004','PF23242'],[]::VARCHAR[],[]::VARCHAR[],1e-5,NULL,NULL,'medium',false,'multi','AAA + TRIP13 lid'),
-- (E) regulation
 ('Ndt80/PhoG family','REG','pfam',['PF05224'],[]::VARCHAR[],[]::VARCHAR[],1e-5,NULL,NULL,'medium',false,'single','family level; Ndt80-like transcription factors are not all meiotic'),
-- family-level counts for the paralog-prone genes (tree in 05-06 separates the meiotic member)
 ('FAMILY kleisin (Rad21/Rec8)','FAMILY','pfam',['PF04824'],[]::VARCHAR[],[]::VARCHAR[],1e-5,NULL,NULL,'high',false,'single','Rad21_Rec8'),
 ('FAMILY MutS (Msh1-6)','FAMILY','pfam',['PF00488'],[]::VARCHAR[],[]::VARCHAR[],1e-5,NULL,NULL,'high',false,'single','MutS_V'),
 ('FAMILY MutL (Mlh1/Pms1/Mlh3)','FAMILY','pfam',['PF08676'],[]::VARCHAR[],[]::VARCHAR[],1e-5,NULL,NULL,'high',false,'single','MutL_C'),
 ('FAMILY HORMA','FAMILY','pfam',['PF02301'],[]::VARCHAR[],[]::VARCHAR[],1e-5,NULL,NULL,'high',false,'single','HORMA (Hop1, Mad2, Rev7)'),
-- ubiquitous CONTROLS, in single-domain and multi-domain style (Ku analysis: multi-domain rules missed 12.5-14% of calls in Kickxellomycota)
 ('CTRL Mre11 (single)','CTRL','pfam',['PF04152'],[]::VARCHAR[],[]::VARCHAR[],1e-5,NULL,NULL,'high',true,'single','Mre11_DNA_bind'),
 ('CTRL Mre11 (multi)','CTRL','pfam',['PF00149','PF04152'],[]::VARCHAR[],[]::VARCHAR[],1e-5,NULL,NULL,'high',true,'multi','nuclease + DNA-binding'),
 ('CTRL Mcm2 (single)','CTRL','pfam',['PF12619'],[]::VARCHAR[],[]::VARCHAR[],1e-5,NULL,NULL,'high',true,'single','MCM2 N-terminal'),
 ('CTRL Mcm2 (multi)','CTRL','pfam',['PF12619','PF00493'],[]::VARCHAR[],[]::VARCHAR[],1e-5,NULL,NULL,'high',true,'multi','MCM2 N-terminal + MCM'),
 ('CTRL Pol delta (single)','CTRL','pfam',['PF03104'],[]::VARCHAR[],[]::VARCHAR[],1e-5,NULL,NULL,'high',true,'single','DNA_pol_B_exo1'),
 ('CTRL Pol delta (multi)','CTRL','pfam',['PF24055','PF00136','PF03104'],[]::VARCHAR[],[]::VARCHAR[],1e-5,NULL,NULL,'high',true,'multi','from the Ku analysis'),
 ('CTRL Pol epsilon (single)','CTRL','pfam',['PF08490'],[]::VARCHAR[],[]::VARCHAR[],1e-5,NULL,NULL,'high',true,'single','DUF1744'),
 ('CTRL Pol epsilon (multi)','CTRL','pfam',['PF08490','PF22634'],[]::VARCHAR[],[]::VARCHAR[],1e-5,NULL,NULL,'high',true,'multi','from the Ku analysis'),
 ('CTRL Rad51 family (single)','CTRL','pfam',['PF08423'],[]::VARCHAR[],[]::VARCHAR[],1e-5,NULL,NULL,'high',true,'single','Rad51 (family level: includes Dmc1 and paralogs)'),
 ('CTRL Rad51 family (multi)','CTRL','pfam',['PF08423','PF14520'],[]::VARCHAR[],[]::VARCHAR[],1e-5,NULL,NULL,'high',true,'multi','from the Ku analysis'),
 ('CTRL MutL (single)','CTRL','pfam',['PF08676'],[]::VARCHAR[],[]::VARCHAR[],1e-5,NULL,NULL,'high',true,'single','mitotic MutL family present everywhere'),
 ('CTRL MutS (single)','CTRL','pfam',['PF00488'],[]::VARCHAR[],[]::VARCHAR[],1e-5,NULL,NULL,'high',true,'single','mitotic MutS family present everywhere'),
 ('CTRL kleisin (single)','CTRL','pfam',['PF04824'],[]::VARCHAR[],[]::VARCHAR[],1e-5,NULL,NULL,'high',true,'single','Rad21 (mitotic cohesin) present everywhere'),
 ('CTRL MutS_I+II (Msh2/3/6)','CTRL','pfam',['PF01624','PF05188'],[]::VARCHAR[],[]::VARCHAR[],1e-5,NULL,NULL,'high',true,'multi','mitotic MutS core; Msh4/Msh5 lack MutS_I'),
 ('CTRL Mlh1 (Mlh1_C single)','CTRL','pfam',['PF16413'],[]::VARCHAR[],[]::VARCHAR[],1e-5,NULL,NULL,'high',true,'single','Mlh1_C, read from the MLH1 references')
) AS t(gene, grp, mode, req_all, req_any, excl, ie_thr, name_re, name_excl_re, conf, is_ctrl, style, note);
