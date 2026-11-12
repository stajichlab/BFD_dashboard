-- Purpose: join genome labels with metrics.
-- Inputs: kx2_h1_genomes, kx2_h1_genome_metrics, kx2_h1_arch_metrics. Output: kx2_h1_genome_all.
-- arch_per1000 = n_arch/ns*1000, where ns = number of Pfam-ANNOTATED proteins in the 2000-protein sample (see 04; ns ranges 53-1867, median 1571),
-- so these are per 1000 annotated proteins, NOT per 1000 sampled proteins. Same for bigram_per1000.
CREATE OR REPLACE TABLE kx2_h1_genome_all AS
SELECT g.*, m.* EXCLUDE(LOCUSTAG), a.* EXCLUDE(LOCUSTAG), a.n_arch::double/a.ns*1000 arch_per1000, a.n_bigram::double/a.ns*1000 bigram_per1000
FROM kx2_h1_genomes g JOIN kx2_h1_genome_metrics m USING(LOCUSTAG) LEFT JOIN kx2_h1_arch_metrics a USING(LOCUSTAG);
