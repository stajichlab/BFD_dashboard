-- 13_coloss.sql
-- Purpose : co-loss summary. For each gene, presence fraction in (A) deep outgroups = Mucoromycota s.l. + root (36 genera),
--           (B) Zoopagales + Entomophthorales + Basidiobolus (12), (C) Ku-RETAINING Kickxellomycota = Dimargaritales +
--           Ramicandelaber (4), (D) the Ku-LOST clade = Kickxellales + Harpellales (backbone node 39; 11 genera),
--           (E) the 107 Coemansia strains; plus the relaxed single-key-domain fraction in (D) where defined (09).
--           Class (parsimony-style reading of the fractions, evaluated on this first-pass tree):
--           'rare everywhere'        : A < 0.25 (gene not reliably detectable even in outgroups; uninformative)
--           'lost with/after Ku'     : A >= 0.5, C >= 0.5, D <= 0.15   (present where Ku is present, gone in the Ku-lost clade)
--           'lost earlier than Ku'   : A >= 0.5, C <= 0.25, B < C or B <= 0.5 (already gone in Ku-retaining Kickxellomycota)
--           'retained'               : A >= 0.5, D >= 0.5
--           'partial / unclear'      : otherwise
-- Inputs  : bfd_scratch.kx2_ku_matrix, kx2_ku_keymatrix, kx2_ku_clade_tips (node 39), kx2_ku_genomes
-- Output  : bfd_scratch.kx2_ku_coloss
-- Thresholds: as in the class definitions above (fractions of genera; strict rule calls from 06).
-- Reports : which genes co-occur with Ku loss, which do not, which were lost earlier.
CREATE OR REPLACE TABLE bfd_scratch.kx2_ku_coloss AS
WITH t AS (
  SELECT m.gene, m.conf, m.category, m.present, m.locustag, m.genus, m.grp, m.set_name,
         CASE WHEN m.set_name='backbone' AND m.grp IN ('Mucoromycota_s.l.','Root') THEN 'A'
              WHEN m.set_name='backbone' AND m.grp IN ('Zoopagales','Entomophthorales_Basidiobolus') THEN 'B'
              WHEN m.set_name='backbone' AND m.grp IN ('Dimargaritales','Ramicandelaber') THEN 'C'
              WHEN m.set_name='backbone' AND m.locustag IN (SELECT tip FROM bfd_scratch.kx2_ku_clade_tips WHERE node=39) THEN 'D'
         END AS blk
  FROM bfd_scratch.kx2_ku_matrix m
),
f AS (
  SELECT gene, conf, category,
    avg(present::INT) FILTER (WHERE blk='A') AS fA, avg(present::INT) FILTER (WHERE blk='B') AS fB,
    avg(present::INT) FILTER (WHERE blk='C') AS fC, avg(present::INT) FILTER (WHERE blk='D') AS fD,
    avg(present::INT) FILTER (WHERE genus='Coemansia') AS fCoem,
    count(*) FILTER (WHERE blk='D') AS nD
  FROM t GROUP BY gene, conf, category
),
k AS (
  SELECT gene, avg(key_present::INT) FILTER (WHERE locustag IN (SELECT tip FROM bfd_scratch.kx2_ku_clade_tips WHERE node=39)) AS keyD,
         avg(key_present::INT) FILTER (WHERE genus='Coemansia') AS keyCoem
  FROM bfd_scratch.kx2_ku_keymatrix GROUP BY gene
)
SELECT f.*, k.keyD, k.keyCoem,
  CASE WHEN fA < 0.25 THEN 'rare everywhere'
       WHEN fA >= 0.5 AND fC >= 0.5 AND fD <= 0.15 THEN 'lost with/after Ku'
       WHEN fA >= 0.5 AND fC <= 0.25 AND (fB <= 0.5 OR fB < fC) THEN 'lost earlier than Ku'
       WHEN fA >= 0.5 AND fD >= 0.5 THEN 'retained'
       ELSE 'partial / unclear' END AS class
FROM f LEFT JOIN k USING (gene);
