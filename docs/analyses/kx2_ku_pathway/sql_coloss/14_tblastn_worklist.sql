-- 14_tblastn_worklist.sql
-- Purpose : worklist for a LATER tblastn/blastx check on the HPCC (NOT run here): genes x genomes where an absence call
--           matters for the Ku question, with SwissProt query accessions of fungal references.
-- Inputs  : bfd_scratch.kx2_ku_matrix, kx2_ku_coloss, kx2_ku_genomes, kx2_ku_dollo, bfd.swissprot_subject (to verify each
--           query accession and carry its gene/protein/organism name)
-- Output  : bfd_scratch.kx2_ku_queries (gene, acc, + names from swissprot_subject)
--           bfd_scratch.kx2_ku_worklist (gene, query accession(s), locustag, genus, ASMID, priority, reason)
-- Selection rules (all must hold):
--   * gene is in the NHEJ / MRN / kinase / telomere list (not the controls) and has coloss class 'lost with/after Ku',
--     'lost earlier than Ku' or 'partial / unclear' with outgroup presence fA >= 0.5 (absence in the ingroup is informative);
--   * the genome is a Kickxellomycota backbone genus or Entomophthorales/Basidiobolus (the lineages around the loss) with
--     present=false for that gene;
--   * plus up to 10 Coemansia strains (highest BUSCO complete, >= 60) with present=false for the gene;
--   * priority 1 = genome inside the Ku-lost clade (node 39) or a Ku-retaining ingroup relative (Dimargaritales, Ramicandelaber)
--     and gene class lost with/after Ku, lost earlier, or Ku itself; priority 2 = everything else on the list.
--   * positive controls (must be found by tblastn): Ku in Ramicandelaber and Dimargaris; Rad51 in all 11 node-39 genera.
CREATE OR REPLACE TABLE bfd_scratch.kx2_ku_queries AS
SELECT q.gene, q.acc, s.gene_name, s.protein_name, s.organism, s.s_len
FROM (VALUES
 ('Ku','P32807'),('Ku','Q04437'),('Ku','Q7SA95'),('Ku','Q7RX73'),('Ku','Q5AVC7'),('Ku','Q5B4H8'),
 ('DNA-PKcs','P78527'),
 ('Lig4','Q5BH83'),('Lig4','Q7SB49'),('Lig4','Q08387'),
 ('XLF/Nej1','Q9P7J2'),('XLF/Nej1','Q06148'),
 ('XRCC4/Lif1','P53150'),('XRCC4/Lif1','Q13426'),
 ('Pol4 (X-family Pol mu/lambda)','Q09693'),('Pol4 (X-family Pol mu/lambda)','P25615'),
 ('Aprataxin (Hnt3)','Q08702'),('Aprataxin (Hnt3)','O74859'),
 ('Sae2/CtIP','P46946'),('Sae2/CtIP','O74986'),
 ('Rif1','P29539'),('Rif1','Q96UP3'),
 ('Taz1/TRF','P79005'),
 ('ATM/Tel1','Q5BHE2'),('ATM/Tel1','O74630'),('ATM/Tel1','P38110'),
 ('ATR/Mec1','Q02099'),('ATR/Mec1','P38111'),
 ('Rad53/Chk2','P22216'),('Rad53/Chk2','Q09170'),
 ('Exo1','P39875'),('Exo1','P53695'),
 ('Nbs1/Xrs2','O43070'),('Nbs1/Xrs2','P33301'),
 ('CTRL Rad51 family','P25454'),('CTRL Rad51 family','Q99133')
) q(gene, acc) LEFT JOIN bfd.swissprot_subject s ON s.accession = q.acc;

CREATE OR REPLACE TABLE bfd_scratch.kx2_ku_worklist AS
WITH genes AS (
  SELECT gene, class, fA FROM bfd_scratch.kx2_ku_coloss
  WHERE gene NOT LIKE 'CTRL%' AND gene NOT LIKE 'WEAK%' AND gene <> 'Ku70' AND gene <> 'Ku80'
    AND (class IN ('lost with/after Ku','lost earlier than Ku') OR (class='partial / unclear' AND fA >= 0.5))
    AND gene IN (SELECT gene FROM bfd_scratch.kx2_ku_queries)
),
ing AS (
  SELECT g.locustag, g.genus, g.ASMID, g.grp, g.set_name, g.busco_complete,
         g.locustag IN (SELECT tip FROM bfd_scratch.kx2_ku_clade_tips WHERE node=39) AS in_ku_lost_clade
  FROM bfd_scratch.kx2_ku_genomes g
  WHERE g.set_name='backbone' AND g.phylum='Kickxellomycota' OR g.grp='Entomophthorales_Basidiobolus'
),
coem AS (
  SELECT locustag, genus, ASMID, grp, set_name, busco_complete, true AS in_ku_lost_clade FROM (
    SELECT g.*, row_number() OVER (ORDER BY busco_complete DESC, n50 DESC) AS rk
    FROM bfd_scratch.kx2_ku_genomes g WHERE g.genus='Coemansia') WHERE rk <= 10
),
tg AS (SELECT * FROM ing UNION ALL SELECT * FROM coem),
cand AS (
  SELECT ge.gene, ge.class, tg.locustag, tg.genus, tg.ASMID, tg.grp, tg.busco_complete, tg.in_ku_lost_clade
  FROM genes ge CROSS JOIN tg
  JOIN bfd_scratch.kx2_ku_matrix m ON m.gene=ge.gene AND m.locustag=tg.locustag AND NOT m.present
)
SELECT c.gene, (SELECT string_agg(q.acc, ',' ORDER BY q.acc) FROM bfd_scratch.kx2_ku_queries q WHERE q.gene=c.gene) AS query_accessions,
       c.locustag, c.genus, c.ASMID, c.grp, c.busco_complete, c.class,
       CASE WHEN (c.in_ku_lost_clade OR c.grp IN ('Dimargaritales','Ramicandelaber')) AND (c.gene='Ku' OR c.class IN ('lost with/after Ku','lost earlier than Ku')) THEN 1 ELSE 2 END AS priority,
       'absent by rule; class=' || c.class AS reason
FROM cand c
UNION ALL
SELECT 'Ku', (SELECT string_agg(acc, ',' ORDER BY acc) FROM bfd_scratch.kx2_ku_queries WHERE gene='Ku'), g.locustag, g.genus, g.ASMID, g.grp, g.busco_complete,
       'terminal single-genus absence', 1, 'Ku absent only in this genus outside the Ku-lost clade (possible annotation miss)'
FROM bfd_scratch.kx2_ku_genomes g
WHERE g.set_name='backbone' AND g.grp NOT IN ('Kickxellales_core','Harpellales','Mycoemilia_Spiromyces','Dimargaritales','Ramicandelaber','Entomophthorales_Basidiobolus')
  AND g.locustag IN (SELECT tip FROM bfd_scratch.kx2_ku_clade_tips WHERE node IN (SELECT node FROM bfd_scratch.kx2_ku_dollo WHERE gene='Ku' AND node <> 39))
UNION ALL
SELECT 'POSITIVE CONTROL Ku', (SELECT string_agg(acc, ',' ORDER BY acc) FROM bfd_scratch.kx2_ku_queries WHERE gene='Ku'), g.locustag, g.genus, g.ASMID, g.grp, g.busco_complete,
       'positive control', 1, 'Ku present by Pfam; tblastn must recover it'
FROM bfd_scratch.kx2_ku_genomes g WHERE g.set_name='backbone' AND g.genus IN ('Ramicandelaber','Dimargaris','Rhizopus')
UNION ALL
SELECT 'POSITIVE CONTROL Rad51', (SELECT string_agg(acc, ',' ORDER BY acc) FROM bfd_scratch.kx2_ku_queries WHERE gene='CTRL Rad51 family'), g.locustag, g.genus, g.ASMID, g.grp, g.busco_complete,
       'positive control', 2, 'Rad51 family present by Pfam in node-39 genera'
FROM bfd_scratch.kx2_ku_genomes g WHERE g.set_name='backbone' AND g.locustag IN (SELECT tip FROM bfd_scratch.kx2_ku_clade_tips WHERE node=39);
