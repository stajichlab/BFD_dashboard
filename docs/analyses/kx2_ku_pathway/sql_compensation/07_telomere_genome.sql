-- Purpose: per-genome telomere features (detection, motif, terminal vs internal, tract length) with assembly-contiguity covariates.
-- Inputs: kx2_kc_wide, bfd.telomere_summary, bfd.telomere_tracts
-- Output: bfd_scratch.kx2_kc_tel_genome (one row per q_loose genome)
-- Thresholds: q_loose from 06; hi_contig = N50 >= 1,000,000 bp ; mid_contig = N50 >= 500,000 bp. Tract length statistics use only tracts of the genome's top monomer.
-- Produces: telomere detection rate per group among contiguous assemblies, motif, median terminal tract length, fraction of tracts that are internal, fraction of telomeric scaffolds with both ends.
-- Caveat: the detector scans canonical G-rich repeats (monomers present in telomere_tracts: CCCTAA, CCCCAA, TTAGGG, CCCTGA, TTGGGG, TCAGGG, ...). A genome with a non-canonical telomere repeat looks like 'no telomere'.
CREATE OR REPLACE TABLE kx2_kc_tel_genome AS
WITH tr AS (
 SELECT LOCUSTAG, count(*) n_tr, median(tract_length) med_len, max(tract_length) max_len,
  median(tract_length) FILTER (terminal) med_len_term, avg((NOT terminal)::int) frac_internal
 FROM bfd.telomere_tracts GROUP BY 1)
SELECT w.LOCUSTAG, w.GENUS, w.ORD, w.CLASS, w.PHYLUM, w.complete_pct, w.N50_bp, w.contig_count, w.n_genes, w.total_length_bp, w.ku,
 coalesce(t.telomere_tracts,0) AS tracts, coalesce(t.telomere_scaffolds,0) AS tel_scaf, coalesce(t.telomere_scaffolds_both_ends,0) AS tel_scaf_both,
 (coalesce(t.telomere_tracts,0)>0) AS has_tel, t.telomere_top_monomer AS monomer, t.telomere_total_length_bp AS tel_bp,
 t.telomere_terminal_count AS n_terminal, t.telomere_internal_count AS n_internal,
 tr.med_len, tr.med_len_term, tr.frac_internal,
 (w.N50_bp>=1000000) AS hi_contig, (w.N50_bp>=500000) AS mid_contig
FROM kx2_kc_wide w LEFT JOIN bfd.telomere_summary t USING (LOCUSTAG) LEFT JOIN tr USING (LOCUSTAG)
WHERE w.q_loose;
SELECT count(*), sum(has_tel::int), sum(hi_contig::int) FROM kx2_kc_tel_genome;
