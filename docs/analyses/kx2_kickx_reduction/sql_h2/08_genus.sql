-- 08_genus.sql
-- Purpose: genus-level medians of per-genome metrics (for Cliff's delta and lineage medians).
-- Inputs: bfd_scratch.kx2_h2_acct, kx2_h2_igorient. Output: bfd_scratch.kx2_h2_genus (filt = 'all' | 'Q').
-- Report numbers: "LINEAGE MEDIANS" table (genes/Mb, coding/intron/intergenic fraction, per-gene bp, masked fraction).
CREATE OR REPLACE TABLE kx2_h2_genus AS
WITH b AS (SELECT a.*, o.cv, o.dv, o.td FROM kx2_h2_acct a LEFT JOIN (SELECT genome, max(med) FILTER (orient='convergent') cv, max(med) FILTER (orient='divergent') dv, max(med) FILTER (orient='tandem') td FROM kx2_h2_igorient GROUP BY 1) o ON o.genome=a.LOCUSTAG WHERE lineage IS NOT NULL),
x AS (SELECT 'all' AS filt,* FROM b UNION ALL SELECT 'Q',* FROM b WHERE q)
SELECT filt, lineage, GENUS, count(*) nsp, median(total_length_bp)/1e6 mb, median(n_genes) genes, median(n_genes/(total_length_bp/1e6)) genes_per_mb,
 median(cds_bp/total_length_bp) coding_frac, median(intron_bp/total_length_bp) intron_frac, median((gap_bp-overlap_bp)/total_length_bp) ig_frac, median(f_intron) f_with_intron, median(mean_span) gene_span,
 median(intron_bp/n_genes) intron_per_gene, median((gap_bp-overlap_bp)/n_genes) ig_per_gene, median(mean_cds) mean_cds, median(masked_bases/total_length_bp) masked_frac, median(ig_med) ig_med, median(f_lt100) f_lt100, median(cv) conv_med, median(dv) div_med, median(td) tand_med, median(contig_count) contigs, median(N50_bp) n50
FROM x GROUP BY ALL;

-- display query for the lineage-medians table:
-- SELECT filt, lineage, count(*) ng, round(median(genes_per_mb)) gpm, round(median(coding_frac),3) cod, round(median(intron_frac),3) intr, round(median(ig_frac),3) ig, round(median(f_with_intron),2) fwi, round(median(gene_span)) span, round(median(intron_per_gene)) ipg, round(median(ig_per_gene)) igpg, round(median(masked_frac),3) mk, round(median(mean_cds)) cds FROM kx2_h2_genus WHERE filt IN ('Q','all') AND lineage IN ('KxCore','Mucorales','Saccharomycotina','Pezizomycotina','Zoopagales','Harpellales','Dimargaritales','Taphrinomycetes','Pneumocystis','Microsporidia','Schizosaccharomyces') GROUP BY ALL ORDER BY filt, lineage;
