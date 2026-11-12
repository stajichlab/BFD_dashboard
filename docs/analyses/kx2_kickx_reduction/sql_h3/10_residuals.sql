-- 10_residuals.sql
-- Purpose : departure of each genus from the kingdom-wide trend: GC3-matched and size-matched residuals,
--           plus within-group scatter on genome size.
-- Inputs  : kx2_h3_genus (all genomes)
-- Outputs : bfd_scratch.kx2_h3_resid (enc_dep, rscu_skew, tai residuals), bfd_scratch.kx2_h3_resid2 (idp, d_outgrp, lowcomp, cys, supply)
-- Method  : reference = non-focal groups (Mucorales, Saccharomycotina other, Pezizomycotina, Basidiomycota, 'other').
--   GC3 bins: floor(gc3*100/5), i.e. 5-percentage-point bins; bin median of reference genera used only if the bin
--   has >=8 reference genera. Size bins: floor(log10(Mb)*6) (~0.17 log10 units), same >=8 rule.
--   residual = genus value - reference bin median. Group value = median of genus residuals.
-- Report  : Kx core ENC-departure residual +0.054, Harpellales -0.017, skew residual +0.095 / -0.029,
--           tai -0.072 / +0.045, gc3 size residual +14.3 / -22.9 points; IDP/outgroup residuals.
CREATE OR REPLACE TABLE bfd_scratch.kx2_h3_resid AS
WITH g AS (SELECT *, floor(gc3*100/5) AS b3, floor(log10(mb)*6) AS bm FROM bfd_scratch.kx2_h3_genus),
ref AS (SELECT * FROM g WHERE hgrp IN ('06 Mucorales','07 Saccharomycotina (other)','14 Pezizomycotina','15 Basidiomycota','99 other')),
t3 AS (SELECT b3, count(*) n3, median(enc_dep) ed3, median(rscu_skew) rs3, median(tai) tai3 FROM ref GROUP BY 1 HAVING count(*)>=8),
tm AS (SELECT bm, count(*) nm, median(enc_dep) edm, median(gc3) gc3m, median(tai) taim FROM ref GROUP BY 1 HAVING count(*)>=8)
SELECT g.hgrp, g.GENUS, g.gc3, g.mb, g.enc_dep, g.enc_dep-t3.ed3 AS enc_dep_res_gc3, g.rscu_skew-t3.rs3 AS skew_res_gc3, g.tai-t3.tai3 AS tai_res_gc3,
  g.enc_dep-tm.edm AS enc_dep_res_size, g.gc3-tm.gc3m AS gc3_res_size, g.tai-tm.taim AS tai_res_size
FROM g LEFT JOIN t3 USING (b3) LEFT JOIN tm USING (bm);

CREATE OR REPLACE TABLE bfd_scratch.kx2_h3_resid2 AS
WITH g AS (SELECT *, floor(gc3*100/5) AS b3 FROM bfd_scratch.kx2_h3_genus),
ref AS (SELECT * FROM g WHERE hgrp IN ('06 Mucorales','07 Saccharomycotina (other)','14 Pezizomycotina','15 Basidiomycota','99 other')),
t3 AS (SELECT b3, median(idp) idp3, median(d_outgrp) do3, median(lowcomp) lc3, median(cys) cys3, median(tai) tai3, median(rscu_supply_r) rs3 FROM ref GROUP BY 1 HAVING count(*)>=8)
SELECT g.hgrp, g.GENUS, g.idp-t3.idp3 idp_res, g.d_outgrp-t3.do3 dout_res, g.lowcomp-t3.lc3 lc_res, g.cys-t3.cys3 cys_res, g.rscu_supply_r-t3.rs3 supply_res
FROM g LEFT JOIN t3 USING (b3);

-- display queries used in the report:
-- SELECT hgrp, count(*) ng, count(enc_dep_res_gc3) n_gc3m, round(median(enc_dep_res_gc3),3) enc_res_gc3, round(median(skew_res_gc3),3) skew_res_gc3,
--   round(median(tai_res_gc3),3) tai_res_gc3, round(median(enc_dep_res_size),3) enc_res_size, round(median(gc3_res_size)*100,1) gc3_res_size_pts,
--   round(median(tai_res_size),3) tai_res_size FROM bfd_scratch.kx2_h3_resid GROUP BY 1 ORDER BY 1;
-- SELECT hgrp, round(median(idp_res)*100,1) idp_res_pts, round(median(dout_res),2) dout_res, round(median(lc_res)*100,2) lowcomp_res_pts,
--   round(median(cys_res)*100,2) cys_res_pts, round(median(supply_res),2) supply_res FROM bfd_scratch.kx2_h3_resid2 GROUP BY 1 ORDER BY 1;
-- size scatter (Pearson r on log10 Mb, genus level; groups with >=4 genera):
-- SELECT hgrp, count(*) n, round(corr(log10(mb), gc3),2) r_gc3_size, round(corr(log10(mb), enc_dep),2) r_encdep_size,
--   round(corr(log10(mb), tai),2) r_tai_size, round(corr(gc3, enc_dep),2) r_gc3_encdep
-- FROM bfd_scratch.kx2_h3_genus GROUP BY ROLLUP(hgrp) HAVING count(*)>=4 ORDER BY 1;
