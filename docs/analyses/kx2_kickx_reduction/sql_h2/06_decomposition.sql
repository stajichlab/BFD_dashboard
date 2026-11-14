-- 06_decomposition.sql
-- Purpose: split each lineage's size deficit vs its baseline into gene-number, CDS/gene, intron/gene, intergenic/gene and remainder terms.
-- Inputs: bfd_scratch.kx2_h2_lin. Result query only (NOT materialised in the original run); output below as a result set.
-- Baselines: KxCore, KxOther, Zoopagales, Harpellales, Dimargaritales vs Mucorales;
--   Saccharomycotina, Taphrinomycetes, Schizosaccharomyces, Pneumocystis, Malassezia, Wallemia, Microsporidia vs Pezizomycotina.
-- Per-gene terms: c = cds, i = intron, g = intergenic, r = remainder, each divided by gene count; pg = c+i+g+r.
-- Ordering A (reported): gene-number effect = (G_base - G_focal) * pg_base; per-gene effects = G_focal * (term_base - term_focal).
-- Ordering B (sensitivity, second query): gene-number effect = (G_base - G_focal) * pg_focal; per-gene effects = G_base * (term_base - term_focal).
-- Report numbers: "ACCOUNTING" table (A) and the "order sensitivity" sentence (B, only filt='Q' and the two key pairs).
-- Q70 rows for Harpellales/Dimargaritales are weak (no Q genomes) and Dimargaritales is not smaller than Mucorales.

-- Ordering A
WITH f AS (SELECT filt, lineage, mb, genes, cds_mb*1e6/genes c, intron_mb*1e6/genes i, ig_mb*1e6/genes g, rem_mb*1e6/genes r FROM kx2_h2_lin),
p AS (SELECT a.filt, a.lineage foc, b.lineage base, b.mb-a.mb defi, (b.genes-a.genes)*(b.c+b.i+b.g+b.r)/1e6 d_genes, a.genes*(b.c-a.c)/1e6 d_cds, a.genes*(b.i-a.i)/1e6 d_intron, a.genes*(b.g-a.g)/1e6 d_ig, a.genes*(b.r-a.r)/1e6 d_rem, b.mb-a.mb - ((b.genes-a.genes)*(b.c+b.i+b.g+b.r)/1e6 + a.genes*((b.c+b.i+b.g+b.r)-(a.c+a.i+a.g+a.r))/1e6) closure
 FROM f a JOIN f b ON a.filt=b.filt AND ((a.lineage IN ('KxCore','KxOther','Zoopagales','Harpellales','Dimargaritales') AND b.lineage='Mucorales') OR (a.lineage IN ('Saccharomycotina','Taphrinomycetes','Schizosaccharomyces','Pneumocystis','Malassezia','Wallemia','Microsporidia') AND b.lineage='Pezizomycotina')))
SELECT filt, foc, base, round(defi,1) defi_mb, round(100*d_genes/defi) pct_genenum, round(100*d_cds/defi) pct_cds_pergene, round(100*d_intron/defi) pct_intron, round(100*d_ig/defi) pct_intergenic, round(100*d_rem/defi) pct_rem, round(closure,2) closure FROM p ORDER BY filt, base, foc;

-- Ordering B (sensitivity; Q filter, the two key pairs)
WITH f AS (SELECT filt, lineage, mb, genes, (cds_mb+intron_mb+ig_mb+rem_mb)*1e6/genes pg, intron_mb*1e6/genes i, ig_mb*1e6/genes g FROM kx2_h2_lin WHERE filt='Q')
SELECT a.lineage, b.lineage base, round(100*(b.genes-a.genes)*a.pg/1e6/(b.mb-a.mb)) pct_genenum_alt, round(100*b.genes*(b.i-a.i)/1e6/(b.mb-a.mb)) pct_intron_alt, round(100*b.genes*(b.g-a.g)/1e6/(b.mb-a.mb)) pct_ig_alt FROM f a JOIN f b ON (a.lineage='KxCore' AND b.lineage='Mucorales') OR (a.lineage='Saccharomycotina' AND b.lineage='Pezizomycotina');
