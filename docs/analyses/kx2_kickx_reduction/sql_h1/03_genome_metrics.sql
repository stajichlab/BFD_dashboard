-- Purpose: per-genome proteome metrics.
-- Inputs: kx2_h1_prot. Output: kx2_h1_genome_metrics.
-- Report: mean/median length, share >1000 aa, fraction >=2 distinct Pfam, any Pfam, domains/100aa, residues in domain,
--         disorder, partial fraction, length-stratified and non-partial multidomain fractions.
CREATE OR REPLACE TABLE kx2_h1_genome_metrics AS
SELECT lt AS LOCUSTAG, count(*) n_prot, sum(length) tot_aa,
 avg(length) mean_len, median(length) med_len,
 avg((length>1000)::int) f_gt1000,
 avg((n_dist>=2)::int) f_multi, avg((n_dom>=1)::int) f_anypfam,
 100.0*sum(n_dom)/sum(length) dom_per100aa,
 sum(dom_res)::double/sum(length) f_res_in_dom,
 avg(is_partial::int) f_partial,
 count(idp_res)::double/count(*) idp_cov,
 sum(idp_res)::double/nullif(sum(length) FILTER (WHERE idp_res IS NOT NULL),0) f_disorder,
 avg((n_dist>=2)::int) FILTER (WHERE length>=300 AND length<600) f_multi_300_600,
 avg((n_dist>=2)::int) FILTER (WHERE length>=600 AND length<1000) f_multi_600_1000,
 avg((n_dist>=2)::int) FILTER (WHERE NOT is_partial) f_multi_nonpartial,
 avg(length) FILTER (WHERE NOT is_partial) mean_len_nonpartial
FROM kx2_h1_prot GROUP BY 1;
