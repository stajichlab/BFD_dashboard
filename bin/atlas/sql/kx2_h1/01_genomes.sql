-- Purpose: define genome set, comparison groups and quality flags for H1 (multidomain/protein-length reduction).
-- Run with default database bfd_scratch (cross-db reads use bfd.<table>).
-- Inputs: bfd_scratch.kw_rep_genomes (7,435 species representatives)
-- Output: bfd_scratch.kx2_h1_genomes
-- Filters: q = N50>=50000 AND (BUSCO complete>=85 OR CLASS Saccharomycetes OR PHYLUM Microsporidia OR GENUS Pneumocystis)
--          q_strict = N50>=50000 AND complete>=85. Rows with NULL GENUS dropped.
-- Report: genome counts per group (core Kickx 149 / 32 q; Mucorales 106 / 78; yeast_other 336 / 311; ...).
CREATE OR REPLACE TABLE kx2_h1_genomes AS
SELECT g.LOCUSTAG, g.ASMID, g.GENUS, g.ORDER_, g.CLASS, g.PHYLUM, g.SPECIES, g.complete_pct, g.N50_bp, g.total_length_bp,
 CASE WHEN g.GENUS IN ('Coemansia','Kickxella','Linderina') THEN 'core_Kickxellales'
      WHEN g.ORDER_='Kickxellales' THEN 'other_Kickxellales'
      WHEN g.ORDER_='Harpellales' THEN 'Harpellales'
      WHEN g.ORDER_='Dimargaritales' THEN 'Dimargaritales'
      WHEN g.ORDER_='Zoopagales' THEN 'Zoopagales'
      WHEN g.ORDER_='Mucorales' THEN 'Mucorales'
      WHEN g.CLASS='Saccharomycetes' AND g.GENUS='Hanseniaspora' THEN 'yeast_Hanseniaspora'
      WHEN g.CLASS='Saccharomycetes' THEN 'yeast_other'
      WHEN g.GENUS='Pneumocystis' THEN 'Pneumocystis'
      WHEN g.ORDER_='Taphrinales' THEN 'Taphrinales'
      WHEN g.PHYLUM='Microsporidia' THEN 'Microsporidia' END AS grp,
 (g.N50_bp>=50000 AND (g.complete_pct>=85 OR g.CLASS='Saccharomycetes' OR g.PHYLUM='Microsporidia' OR g.GENUS='Pneumocystis')) AS q,
 (g.N50_bp>=50000 AND g.complete_pct>=85) AS q_strict
FROM kw_rep_genomes g WHERE g.GENUS IS NOT NULL AND (g.ORDER_ IN ('Kickxellales','Harpellales','Dimargaritales','Zoopagales','Mucorales','Taphrinales') OR g.CLASS='Saccharomycetes' OR g.GENUS='Pneumocystis' OR g.PHYLUM='Microsporidia');
-- check query
SELECT grp, count(*) n, count(distinct GENUS) g, sum(q::int) nq, sum(q_strict::int) nqs FROM kx2_h1_genomes GROUP BY 1 ORDER BY 1;
