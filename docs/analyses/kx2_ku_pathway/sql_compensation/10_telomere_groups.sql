-- Purpose: telomere detection, motif, terminal tract length per clade group and per Ku status, restricted to contiguous assemblies; length rank of Harpellales genera.
-- Inputs: kx2_kc_tel_genome (07), kx2_kc_grp (08)
-- Outputs: bfd_scratch.kx2_kc_tel_group (group x Ku-status), bfd_scratch.kx2_kc_tel_genus (genus level, mid-contig genomes), bfd_scratch.kx2_kc_tel_lenrank (terminal tract length percentile among genera)
-- Thresholds: mid_contig = N50 >= 500 kb (genus-level detection and length); a genus is 'detected' if >= 50% of its mid-contig genomes have >= 1 tract.
-- Produces: telomere detection rate among contiguous assemblies (Kickxellales 0 of 1 genus (Coemansia, 2 genomes); Harpellales 3 of 3 genera; Mucoromycota+Mortierellomycota+Glomeromycota ...), motif CCCTAA, terminal tract length, Zancudomyces rank.
-- Note: telomere_summary reports terminal tracts only (telomere_internal_count = 0 for all 23,582 genomes), so interstitial repeats cannot be assessed from this table.
CREATE OR REPLACE TABLE kx2_kc_tel_genus AS
SELECT t.PHYLUM, t.GENUS, g.grp, count(*) AS n_mid, avg(t.has_tel::int) AS det, median(t.tracts) FILTER (t.has_tel) AS tracts,
 median(t.med_len_term) AS len_term, median(t.tel_scaf_both::double/nullif(t.tel_scaf,0)) AS both_frac, mode(t.monomer) AS monomer, median(t.N50_bp) AS n50, any_value(t.ku) AS ku_example
FROM kx2_kc_tel_genome t JOIN kx2_kc_grp g USING (PHYLUM, GENUS) WHERE t.mid_contig GROUP BY ALL;

CREATE OR REPLACE TABLE kx2_kc_tel_group AS
SELECT grp, count(*) AS genera_midcontig, sum((det>=0.5)::int) AS genera_detected, avg(det) AS mean_det_rate, median(len_term) AS med_terminal_len, median(tracts) AS med_tracts, median(both_frac) AS med_both_frac, median(n50) AS med_n50
FROM kx2_kc_tel_genus GROUP BY 1;

CREATE OR REPLACE TABLE kx2_kc_tel_lenrank AS
WITH g AS (SELECT PHYLUM, GENUS, median(med_len_term) AS len FROM kx2_kc_tel_genome WHERE has_tel AND mid_contig GROUP BY 1,2)
SELECT g.GENUS, g.len, (SELECT avg((g2.len<=g.len)::int) FROM g g2) AS pctile, (SELECT count(*) FROM g) AS n_genera
FROM g WHERE g.GENUS IN ('Zancudomyces','Smittium','Capniomyces');
SELECT * FROM kx2_kc_tel_group ORDER BY 1;
SELECT * FROM kx2_kc_tel_lenrank;
