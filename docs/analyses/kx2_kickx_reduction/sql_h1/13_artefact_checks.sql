-- Purpose: artefact controls (result queries only, no tables).
-- Inputs: kx2_h1_genome_metrics, kx2_h1_genomes, kx2_h1_prot.
-- (a) genome-level correlation of multidomain fraction / mean length with partial fraction, N50, BUSCO (ALL genomes of the group, not only q).
-- (b) non-partial proteins in q genomes stratified by length bin: multidomain and any-Pfam fraction (pooled proteins).
-- Report: artefact checks section (r = -0.70 multi vs partial in core; length-stratified table).
WITH x AS (SELECT g.grp, g.GENUS, f.LOCUSTAG, f.n_prot, f.f_multi, f.f_partial, g.N50_bp, g.complete_pct, f.mean_len FROM kx2_h1_genome_metrics f JOIN kx2_h1_genomes g USING(LOCUSTAG) WHERE g.grp IN ('core_Kickxellales','Mucorales','yeast_other'))
SELECT grp, count(*) n, round(corr(f_multi, f_partial),2) r_multi_partial, round(corr(f_multi, ln(N50_bp)),2) r_multi_lnN50, round(corr(f_multi, complete_pct),2) r_multi_busco, round(corr(mean_len, f_partial),2) r_len_partial, round(corr(mean_len, ln(N50_bp)),2) r_len_N50, round(corr(f_multi, mean_len),2) r_multi_len FROM x GROUP BY 1;

SELECT grp, CASE WHEN length<300 THEN '1_<300' WHEN length<600 THEN '2_300-600' WHEN length<1000 THEN '3_600-1000' WHEN length<2000 THEN '4_1000-2000' ELSE '5_>=2000' END lb, count(*) n, round(avg((n_dist>=2)::int),3) f_multi, round(avg((n_dom>=1)::int),3) f_any, round(100.0*count(*)/sum(count(*)) OVER (PARTITION BY grp),2) pct_prots FROM kx2_h1_prot p JOIN kx2_h1_genomes g ON g.LOCUSTAG=p.lt WHERE g.q AND g.grp IN ('core_Kickxellales','Mucorales','yeast_other') AND NOT p.is_partial GROUP BY 1,2 ORDER BY 2,1;
