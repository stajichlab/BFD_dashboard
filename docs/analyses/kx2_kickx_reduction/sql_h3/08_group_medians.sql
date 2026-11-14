-- 08_group_medians.sql
-- Purpose : per-group medians of genus medians (filt='all' = all reps; filt='b85' = BUSCO>=85 reps).
-- Inputs  : kx2_h3_genus, kx2_h3_genus_f
-- Output  : bfd_scratch.kx2_h3_grpmed
-- Notes   : gc3, fymink, asn, lys, cys, gln, idp are x100 (percent). Microsporidia, Harpellales,
--           Schizosaccharomyces and fast-evolving Sacch have no b85 rows.
-- Report  : the per-group medians table (e.g. Kx core: gc3 67.2, enc 49.5, enc_dep 0.110, tai 0.287).
CREATE OR REPLACE TABLE bfd_scratch.kx2_h3_grpmed AS
SELECT 'all' AS filt, hgrp, count(*) ng, sum(n_sp) nsp, round(median(mb),1) mb, round(median(gc),1) gc, round(median(gc3)*100,1) gc3, round(median(enc),1) enc, round(median(enc_dep),3) enc_dep, round(median(rscu_skew),3) rscu_skew,
 median(n_trna) n_trna, median(n_antic) n_antic, round(median(tai),3) tai, round(median(rscu_supply_r),2) r_supply, round(median(fymink)*100,1) fymink, round(median(asn)*100,2) asn, round(median(lys)*100,2) lys, round(median(cys)*100,2) cys, round(median(gln)*100,2) gln, round(median(idp)*100,1) idp,
 round(median(d_schpo),2) d_schpo, round(median(f60_schpo),3) f60, round(median(cov_schpo),2) cov, round(median(d_outgrp),2) d_out
FROM bfd_scratch.kx2_h3_genus GROUP BY 2
UNION ALL
SELECT 'b85', hgrp, count(*), sum(n_sp), round(median(mb),1), round(median(gc),1), round(median(gc3)*100,1), round(median(enc),1), round(median(enc_dep),3), round(median(rscu_skew),3),
 median(n_trna), median(n_antic), round(median(tai),3), round(median(rscu_supply_r),2), round(median(fymink)*100,1), round(median(asn)*100,2), round(median(lys)*100,2), round(median(cys)*100,2), round(median(gln)*100,2), round(median(idp)*100,1),
 round(median(d_schpo),2), round(median(f60_schpo),3), round(median(cov_schpo),2), round(median(d_outgrp),2)
FROM bfd_scratch.kx2_h3_genus_f GROUP BY 2;

-- display query (per-group table in the report):
-- SELECT * FROM bfd_scratch.kx2_h3_grpmed WHERE filt='all' ORDER BY hgrp;
