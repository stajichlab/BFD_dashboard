-- 05_lin.sql
-- Purpose: lineage-level accounting terms (Mb). Genus median of each component, then mean across genera of a lineage.
-- Inputs: bfd_scratch.kx2_h2_acct.
-- Output: bfd_scratch.kx2_h2_lin  (filt = 'all' | 'Q' | 'Q70')
-- Filters: Q = busco_complete>=85 AND N50_bp>=50000; Q70 = busco_complete>=70 AND N50_bp>=25000; 'all' = no filter.
--   ig_net = positive gaps minus overlaps. rem = remainder_bp.
-- Report numbers: "Lineage medians" table and the Mb inputs for the decomposition (06).
CREATE OR REPLACE TABLE kx2_h2_lin AS
WITH base AS (SELECT *, (gap_bp-overlap_bp) ig_net, 'all' AS filt FROM kx2_h2_acct WHERE lineage IS NOT NULL
 UNION ALL SELECT *, (gap_bp-overlap_bp), 'Q' FROM kx2_h2_acct WHERE lineage IS NOT NULL AND q
 UNION ALL SELECT *, (gap_bp-overlap_bp), 'Q70' FROM kx2_h2_acct WHERE lineage IS NOT NULL AND busco_complete>=70 AND N50_bp>=25000),
gm AS (SELECT filt, lineage, GENUS, count(*) nsp, median(total_length_bp) tot, median(n_genes) genes, median(cds_bp) cds, median(intron_bp) intr, median(ig_net) ig, median(remainder_bp) rem,
  median(masked_bases) masked, median(total_n_bases) nb, median(f_intron) f_intron, median(ipt) ipt, median(mean_cds) mean_cds, median(mean_span) mean_span, median(ig_med) ig_med, median(f_lt100) f_lt100, median(f_ov) f_ov, median(contig_count) contigs
  FROM base GROUP BY ALL)
SELECT filt, lineage, count(*) n_genera, sum(nsp) n_sp, avg(tot)/1e6 mb, avg(genes) genes, avg(cds)/1e6 cds_mb, avg(intr)/1e6 intron_mb, avg(ig)/1e6 ig_mb, avg(rem)/1e6 rem_mb, avg(masked)/1e6 masked_mb, avg(nb)/1e6 n_mb,
 avg(f_intron) f_intron, avg(ipt) ipt, avg(mean_cds) mean_cds, avg(mean_span) mean_span, avg(ig_med) ig_med, avg(f_lt100) f_lt100, avg(f_ov) f_ov, avg(contigs) contigs
FROM gm GROUP BY 1,2;

-- display query used for the report (closure = mb - cds - intron - ig - rem):
-- SELECT filt, lineage, n_genera, n_sp, round(mb,1) mb, round(genes) genes, round(cds_mb,2) cds, round(intron_mb,2) intr, round(ig_mb,2) ig, round(rem_mb,2) rem, round(mb-cds_mb-intron_mb-ig_mb-rem_mb,3) closure, round(masked_mb,2) masked, round(f_intron,2) fi, round(ig_med) igmed, round(f_lt100,3) lt100, round(f_ov,4) fov FROM kx2_h2_lin ORDER BY filt, lineage;
