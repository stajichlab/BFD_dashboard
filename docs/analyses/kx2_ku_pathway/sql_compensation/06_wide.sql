-- Purpose: zero-filled genome x feature matrix with quality flags; feature map documents how genes roll up into features.
-- Inputs: kx2_kc_genomes, kx2_kc_gene_counts
-- Outputs: bfd_scratch.kx2_kc_featmap (gene -> feature), bfd_scratch.kx2_kc_wide (one row per genome with BUSCO, N50, n_genes, flags and feature counts)
-- Thresholds: q_loose = BUSCO complete >= 60 AND n_genes >= 3000 AND BUSCO duplicated <= 10 (same filter for every group, because early-diverging
--             fungi have BUSCO ~70-85%); q_strict = BUSCO complete >= 85 AND duplicated <= 10 (sensitivity). Genomes with NULL GENUS/BUSCO dropped.
-- Produces: the unit table for all group/genus comparisons.
CREATE OR REPLACE TABLE kx2_kc_featmap AS SELECT * FROM (VALUES
 ('KU_pf','ku'),('POLQ_pf','polq'),('PARP1like_pf','parp1l'),('PARP_any_pf','parp_any'),
 ('LIG1','lig1'),('LIG3','lig3'),('LIG4','lig4'),('fam_ligA_Mdom','ligfam'),
 ('FEN1','fen1'),('EXO1','exo1'),('GEN1','gen1'),('fam_XPGdom','xpgfam'),
 ('RAD51','rad51'),('DMC1','dmc1'),
 ('RAD51_PARALOG_BCD','rad51para'),('RAD55','rad51para'),('RAD57','rad51para'),('XRCC2','rad51para'),('XRCC3','rad51para'),
 ('fam_Rad51dom','rad51fam'),('RAD52','rad52_59'),('RAD59','rad52_59'),('fam_Rad52dom','rad52fam'),
 ('fam_Rad54Ndom','rad54N'),('RAD54','rad54name'),('RDH54','rad54name'),('RAD54B','rad54name'),
 ('BRCA2','brca2name'),('fam_BRCA2dom','brca2dom'),
 ('MUS81','mus81'),('MMS4_EME1','mms4'),('XPF_RAD1','xpf'),('RAD10_ERCC1','rad10name'),('fam_Rad10dom','rad10dom'),('fam_ERCC4dom','ercc4fam'),
 ('TERT','tert'),('fam_TERTdom','tertdom'),('XRCC1','xrcc1'),('fam_DNApolAdom','polAfam')
) t(gene, feature);

CREATE OR REPLACE TABLE kx2_kc_wide AS
WITH g AS (SELECT * FROM kx2_kc_genomes WHERE GENUS IS NOT NULL AND complete_pct IS NOT NULL),
f AS (SELECT DISTINCT feature FROM kx2_kc_featmap),
agg AS (SELECT c.LOCUSTAG, m.feature, sum(c.n) AS n FROM kx2_kc_gene_counts c JOIN kx2_kc_featmap m USING (gene) GROUP BY 1,2),
long AS (SELECT g.LOCUSTAG, f.feature, coalesce(a.n,0) AS n FROM g CROSS JOIN f LEFT JOIN agg a ON a.LOCUSTAG=g.LOCUSTAG AND a.feature=f.feature),
p AS (PIVOT long ON feature USING sum(n) GROUP BY LOCUSTAG)
SELECT g.*, (g.complete_pct>=60 AND g.n_genes>=3000 AND g.duplicated_pct<=10) AS q_loose, (g.complete_pct>=85 AND g.duplicated_pct<=10) AS q_strict, p.* EXCLUDE (LOCUSTAG)
FROM g JOIN p USING (LOCUSTAG);
SELECT count(*), sum(q_loose::int), sum(q_strict::int), sum((q_loose AND ku>0)::int), sum((q_loose AND ku=0)::int) FROM kx2_kc_wide;
