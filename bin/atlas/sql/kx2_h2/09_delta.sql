-- 09_delta.sql
-- Purpose: Cliff's delta between genera of focal and reference lineages, for 17 metrics.
-- Inputs: bfd_scratch.kx2_h2_genus. Output: bfd_scratch.kx2_h2_delta (+ pivot result query).
-- Pairs: KxCore vs Mucorales, KxCore vs Saccharomycotina, Saccharomycotina vs Pezizomycotina, Harpellales vs Mucorales.
--   Cliff's delta = mean(sign(focal - ref)) over all genus pairs. Harpellales has no Q genomes, so only filt='all' is meaningful.
-- Report numbers: "Cliff's delta (genus level, Q)" list.
CREATE OR REPLACE TABLE kx2_h2_delta AS
WITH m AS (SELECT filt, lineage, GENUS, metric, v FROM kx2_h2_genus UNPIVOT (v FOR metric IN (mb, genes, genes_per_mb, coding_frac, intron_frac, ig_frac, f_with_intron, gene_span, intron_per_gene, ig_per_gene, mean_cds, masked_frac, ig_med, f_lt100, conv_med, div_med, tand_med)))
SELECT a.filt, a.lineage AS foc, b.lineage AS reflin, a.metric, count(distinct a.GENUS) nf, count(distinct b.GENUS) nr, median(a.v) med_f, round(avg(sign(a.v-b.v)),2) cliffs_delta
FROM m a JOIN m b ON a.filt=b.filt AND a.metric=b.metric AND ((a.lineage='KxCore' AND b.lineage IN ('Mucorales','Saccharomycotina')) OR (a.lineage='Saccharomycotina' AND b.lineage='Pezizomycotina') OR (a.lineage='Harpellales' AND b.lineage='Mucorales'))
GROUP BY ALL;

-- pivot result query:
-- SELECT metric, max(CASE WHEN foc='KxCore' AND reflin='Mucorales' THEN cliffs_delta END) kx_vs_muc, max(CASE WHEN foc='Saccharomycotina' AND reflin='Pezizomycotina' THEN cliffs_delta END) sac_vs_pez, max(CASE WHEN foc='KxCore' AND reflin='Saccharomycotina' THEN cliffs_delta END) kx_vs_sac FROM kx2_h2_delta WHERE filt='Q' GROUP BY 1 ORDER BY 1;
