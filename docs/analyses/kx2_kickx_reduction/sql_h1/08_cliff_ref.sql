-- Purpose: Cliff's delta core Kickxellales (filtered) vs unreduced Zoopagomycota/Kickxellomycota genera
--          (unfiltered: Harpellales, Dimargaritales, Zoopagales, other Kickxellales; 19 genera).
-- Inputs: kx2_h1_genus. Output: kx2_h1_cliff_ref. Report: "core vs the unreduced ... (19 genera)" deltas.
CREATE OR REPLACE TABLE kx2_h1_cliff_ref AS
WITH L AS (SELECT mode, grp, GENUS, m, v FROM kx2_h1_genus UNPIVOT(v FOR m IN (mean_len, f_gt1000, f_multi, f_anypfam, dom_per100aa, f_res_in_dom, arch_per1000, bigram_per1000, promisc, f_multi_np)) WHERE v IS NOT NULL),
A AS (SELECT * FROM L WHERE mode='filtered' AND grp='core_Kickxellales'),
R AS (SELECT * FROM L WHERE mode='unfiltered' AND grp IN ('Harpellales','Dimargaritales','Zoopagales','other_Kickxellales'))
SELECT A.m, count(distinct R.GENUS) n_ref_genera, round(sum(sign(A.v-R.v))/count(*),2) cliff_core_vs_ref, round(median(R.v),3) ref_median, round(median(A.v),3) core_median FROM A JOIN R USING(m) GROUP BY 1;
SELECT * FROM kx2_h1_cliff_ref ORDER BY 1;
