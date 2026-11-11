-- Purpose: Cliff's delta of each group vs Mucorales over genus medians.
-- Inputs: kx2_h1_genus. Output: kx2_h1_cliff.
-- delta = (#pairs a>b - #pairs a<b)/#pairs; genera are the units. Report: delta table (filtered mode).
CREATE OR REPLACE TABLE kx2_h1_cliff AS
WITH L AS (SELECT mode, grp, GENUS, m, v FROM kx2_h1_genus UNPIVOT(v FOR m IN (mean_len, med_len, f_gt1000, f_multi, f_anypfam, dom_per100aa, f_res_in_dom, f_disorder, arch_per1000, bigram_per1000, promisc, f_multi_300_600, f_multi_600_1000, f_multi_np, mean_len_np)) WHERE v IS NOT NULL),
R AS (SELECT * FROM L WHERE grp='Mucorales')
SELECT a.mode, a.grp, a.m, count(distinct a.GENUS) n_gen, round((sum(sign(a.v-r.v)))/ (count(*)),2) cliff_vs_Mucorales
FROM L a JOIN R r ON a.mode=r.mode AND a.m=r.m WHERE a.grp<>'Mucorales' GROUP BY ALL;
SELECT m, max(cliff_vs_Mucorales) FILTER (WHERE grp='core_Kickxellales') coreK, max(cliff_vs_Mucorales) FILTER (WHERE grp='Harpellales') Harp, max(cliff_vs_Mucorales) FILTER (WHERE grp='Dimargaritales') Dim, max(cliff_vs_Mucorales) FILTER (WHERE grp='yeast_other') yeast, max(cliff_vs_Mucorales) FILTER (WHERE grp='Microsporidia') Micro, max(cliff_vs_Mucorales) FILTER (WHERE grp='Taphrinales') Taph FROM kx2_h1_cliff WHERE mode='filtered' GROUP BY 1 ORDER BY 1;
