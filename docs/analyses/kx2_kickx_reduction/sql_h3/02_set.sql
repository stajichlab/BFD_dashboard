-- 02_set.sql
-- Purpose : genome set and H3 group assignment (species representatives with codon data).
-- Inputs  : bfd_scratch.kx_genomes (is_sp_rep; 1,919 genomes), bfd_scratch.kw_rep_genomes,
--           bfd.asm_stats (gc_pct, masked_pct for the added yeasts), bfd.codon_freq (presence filter)
-- Output  : bfd_scratch.kx2_h3_set  (1,926 genomes)
-- Notes   : Hanseniaspora/Saccharomycodes are replaced by all their reps from kw_rep_genomes
--           (kx_genomes had 1 Hanseniaspora; kw_rep_genomes has 23+2). Group label hgrp is the
--           report group. 'Kx core' = Kx:Coemansia + Kx:Kickxellales-core-other.
--           BUSCO>=85 variant is applied later (busco_complete>=85) in 08/09.
-- Report  : group n / genus counts in the table; BUSCO>=85 counts (0 for Harpellales, Microsporidia,
--           Schizosaccharomyces, fast-evolving Sacch).
CREATE OR REPLACE TABLE bfd_scratch.kx2_h3_set AS
SELECT LOCUSTAG AS g, GENUS, grp, total_length_bp/1e6 AS mb, gc_pct, busco_complete, N50_bp, masked_pct,
 CASE
  WHEN grp IN ('Kx:Coemansia','Kx:Kickxellales-core-other') THEN '01 Kx core Kickxellales'
  WHEN grp='Kx:Kickxellales-early' THEN '02 Kx early Kickxellales'
  WHEN grp='Kx:Harpellales' THEN '03 Harpellales'
  WHEN grp='Kx:Dimargaritales' THEN '04 Dimargaritales'
  WHEN grp='Zoopagales' THEN '05 Zoopagales'
  WHEN grp='Mucorales' THEN '06 Mucorales'
  WHEN GENUS IN ('Hanseniaspora','Saccharomycodes') THEN '08 Sacch fast-evolving (Hanseniaspora/Saccharomycodes)'
  WHEN grp LIKE 'Sacch:%' THEN '07 Saccharomycotina (other)'
  WHEN grp='Microsporidia' THEN '09 Microsporidia'
  WHEN grp='Pneumocystis' THEN '10 Pneumocystis'
  WHEN grp='Taphrinomycetes' THEN '11 Taphrinomycetes'
  WHEN grp='Schizosaccharomyces' THEN '12 Schizosaccharomyces'
  WHEN grp='Malassezia' THEN '13 Malassezia'
  WHEN grp LIKE 'Pez:%' THEN '14 Pezizomycotina'
  WHEN grp IN ('Agaricomycetes','Exobasidiomycetes','Microbotryomycetes','Pucciniomycetes','Tremellomycetes','Ustilaginomycetes','Wallemia') THEN '15 Basidiomycota'
  ELSE '99 other' END AS hgrp
FROM bfd_scratch.kx_genomes WHERE is_sp_rep AND LOCUSTAG IN (SELECT distinct species_prefix FROM bfd.codon_freq);

-- replace the (single) kx_genomes Hanseniaspora/Saccharomycodes rows by all kw_rep_genomes reps
DELETE FROM bfd_scratch.kx2_h3_set WHERE hgrp LIKE '08%';
INSERT INTO bfd_scratch.kx2_h3_set
SELECT r.LOCUSTAG, r.GENUS, 'Sacch:Saccharomycetes', r.total_length_bp/1e6, a.gc_pct, r.complete_pct, r.N50_bp, a.masked_pct,
 '08 Sacch fast-evolving (Hanseniaspora/Saccharomycodes)'
FROM bfd_scratch.kw_rep_genomes r LEFT JOIN bfd.asm_stats a ON a.ASMID=r.ASMID
WHERE r.GENUS IN ('Hanseniaspora','Saccharomycodes') AND r.LOCUSTAG IN (SELECT distinct species_prefix FROM bfd.codon_freq);
