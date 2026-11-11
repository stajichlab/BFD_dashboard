-- 09_cliff.sql
-- Purpose : Cliff's delta (genus level) for focal H3 groups vs baseline groups, both filters.
-- Inputs  : kx2_h3_genus, kx2_h3_genus_f
-- Outputs : bfd_scratch.kx2_h3_long (long metric table), bfd_scratch.kx2_h3_cliff
-- Definition: delta = (#(x>y) - #(x<y)) / (nx*ny) over all genus pairs (x focal, y reference); sign(x-y)
--   summed over the cross join. +1 means every focal genus exceeds every reference genus.
-- Focal   : Kx core, Harpellales, Dimargaritales, Zoopagales, fast-evolving Sacch, Microsporidia, Pneumocystis.
-- Refs    : Mucorales, Saccharomycotina (other), Pezizomycotina, Basidiomycota.
-- Caution : Kx core is 5 genera (4 with BUSCO>=85); Harpellales 4 (0 with BUSCO>=85). No p-values.
-- Report  : e.g. Kx core vs Mucorales enc_dep +0.87, gc3 +0.97, tai -0.72, outgroup identity -1.
CREATE OR REPLACE TABLE bfd_scratch.kx2_h3_long AS
WITH u AS (
 SELECT 'all' filt, * FROM bfd_scratch.kx2_h3_genus
 UNION ALL BY NAME SELECT 'b85' filt, * FROM bfd_scratch.kx2_h3_genus_f)
SELECT filt, hgrp, GENUS, metric, val FROM (
 UNPIVOT (SELECT filt, hgrp, GENUS, mb, gc, gc3, enc, enc_dep, rscu_skew, n_trna, trna_per_mb, tai, rscu_supply_r, fymink, asn, lys, cys, gln, lowcomp, idp, d_outgrp, d_schpo FROM u)
 ON mb, gc, gc3, enc, enc_dep, rscu_skew, n_trna, trna_per_mb, tai, rscu_supply_r, fymink, asn, lys, cys, gln, lowcomp, idp, d_outgrp, d_schpo INTO NAME metric VALUE val)
WHERE val IS NOT NULL;

CREATE OR REPLACE TABLE bfd_scratch.kx2_h3_cliff AS
SELECT x.filt, x.hgrp AS focal, y.hgrp AS ref, x.metric, count(distinct x.GENUS) nx, count(distinct y.GENUS) ny,
  round(sum(sign(x.val-y.val))/ (count(distinct x.GENUS)*count(distinct y.GENUS)),2) AS cliffs_delta
FROM bfd_scratch.kx2_h3_long x JOIN bfd_scratch.kx2_h3_long y ON x.filt=y.filt AND x.metric=y.metric
WHERE x.hgrp IN ('01 Kx core Kickxellales','03 Harpellales','04 Dimargaritales','05 Zoopagales','08 Sacch fast-evolving (Hanseniaspora/Saccharomycodes)','09 Microsporidia','10 Pneumocystis')
 AND y.hgrp IN ('06 Mucorales','07 Saccharomycotina (other)','14 Pezizomycotina','15 Basidiomycota')
GROUP BY ALL;

-- display query used in the report:
-- SELECT metric, filt, focal, ref, nx, ny, cliffs_delta FROM bfd_scratch.kx2_h3_cliff
-- WHERE focal IN ('01 Kx core Kickxellales','03 Harpellales')
--   AND metric IN ('enc','enc_dep','gc3','rscu_skew','tai','n_trna','fymink','idp','d_outgrp','rscu_supply_r')
--   AND ref IN ('06 Mucorales','07 Saccharomycotina (other)') ORDER BY metric, filt, focal, ref;
