-- Purpose: genus medians (replication unit), filtered and unfiltered.
-- Inputs: kx2_h1_genome_all. Output: kx2_h1_genus.
-- 'filtered' = q genomes plus ALL Harpellales/Dimargaritales genomes (none pass q); 'unfiltered' = all genomes.
-- Report: the per-group medians table (second query = median over genera per group).
CREATE OR REPLACE TABLE kx2_h1_genus AS
WITH sel AS (SELECT 'filtered' AS mode, * FROM kx2_h1_genome_all WHERE q OR grp IN ('Harpellales','Dimargaritales')
 UNION ALL SELECT 'unfiltered', * FROM kx2_h1_genome_all)
SELECT mode, grp, GENUS, count(*) n_genomes, median(n_prot) n_prot, median(mean_len) mean_len, median(med_len) med_len, median(f_gt1000) f_gt1000, median(f_multi) f_multi, median(f_anypfam) f_anypfam, median(dom_per100aa) dom_per100aa, median(f_res_in_dom) f_res_in_dom, median(f_disorder) f_disorder, median(f_partial) f_partial, median(arch_per1000) arch_per1000, median(bigram_per1000) bigram_per1000, median(promisc) promisc, median(f_multi_300_600) f_multi_300_600, median(f_multi_600_1000) f_multi_600_1000, median(f_multi_nonpartial) f_multi_np, median(mean_len_nonpartial) mean_len_np, median(complete_pct) busco, median(N50_bp) n50, median(idp_cov) idp_cov
FROM sel GROUP BY ALL;
SELECT mode, grp, count(*) ngen, round(median(n_prot)) nprot, round(median(mean_len)) mean_len, round(median(med_len)) med_len, round(median(f_gt1000),3) gt1000, round(median(f_multi),3) multi, round(median(f_anypfam),3) anyp, round(median(dom_per100aa),2) d100, round(median(f_res_in_dom),3) resdom, round(median(f_disorder),3) disord, round(median(arch_per1000),0) arch1000, round(median(bigram_per1000),0) bg1000, round(median(promisc),2) promisc, round(median(f_partial),3) part FROM kx2_h1_genus GROUP BY 1,2 ORDER BY 1,2;
