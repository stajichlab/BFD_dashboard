-- Purpose: long table of gene/family counts per genome (distinct genes: -T<n> isoform suffix stripped).
-- Inputs: kx2_kc_calls
-- Output: bfd_scratch.kx2_kc_gene_counts (LOCUSTAG, gene, n); absent genome x gene = no row (zero-filled in 06).
-- Rules: gene = class (cls) when pf_ok; Pfam-defined sets: KU_pf, POLQ_pf, PARP1like_pf, PARP_any_pf; family sets fam_* count every protein carrying the domain
--        (independent of best-hit naming; the primary evidence for divergent families).
-- Produces: all copy-number / presence numbers in the report.
CREATE OR REPLACE TABLE kx2_kc_gene_counts AS
WITH u AS (
 SELECT LOCUSTAG, protein_id, cls AS gene FROM kx2_kc_calls WHERE cls IS NOT NULL AND pf_ok
 UNION ALL SELECT LOCUSTAG, protein_id, 'KU_pf' FROM kx2_kc_calls WHERE ku_pf
 UNION ALL SELECT LOCUSTAG, protein_id, 'POLQ_pf' FROM kx2_kc_calls WHERE polq_pf
 UNION ALL SELECT LOCUSTAG, protein_id, 'PARP_any_pf' FROM kx2_kc_calls WHERE parp_pf
 UNION ALL SELECT LOCUSTAG, protein_id, 'PARP1like_pf' FROM kx2_kc_calls WHERE parp1like_pf
 UNION ALL SELECT LOCUSTAG, protein_id, 'fam_Rad51dom' FROM kx2_kc_calls WHERE list_contains(doms,'Rad51')
 UNION ALL SELECT LOCUSTAG, protein_id, 'fam_Rad52dom' FROM kx2_kc_calls WHERE list_contains(doms,'Rad52_Rad22')
 UNION ALL SELECT LOCUSTAG, protein_id, 'fam_Rad54Ndom' FROM kx2_kc_calls WHERE list_contains(doms,'Rad54_N')
 UNION ALL SELECT LOCUSTAG, protein_id, 'fam_Rad10dom' FROM kx2_kc_calls WHERE list_contains(doms,'Rad10')
 UNION ALL SELECT LOCUSTAG, protein_id, 'fam_BRCA2dom' FROM kx2_kc_calls WHERE list_contains(doms,'BRCA2') OR list_contains(doms,'BRCA-2_OB1')
 UNION ALL SELECT LOCUSTAG, protein_id, 'fam_ERCC4dom' FROM kx2_kc_calls WHERE list_contains(doms,'ERCC4')
 UNION ALL SELECT LOCUSTAG, protein_id, 'fam_XPGdom' FROM kx2_kc_calls WHERE list_contains(doms,'XPG_I')
 UNION ALL SELECT LOCUSTAG, protein_id, 'fam_ligA_Mdom' FROM kx2_kc_calls WHERE list_contains(doms,'DNA_ligase_A_M')
 UNION ALL SELECT LOCUSTAG, protein_id, 'fam_TERTdom' FROM kx2_kc_calls WHERE list_contains(doms,'Telomerase_RBD')
 UNION ALL SELECT LOCUSTAG, protein_id, 'fam_DNApolAdom' FROM kx2_kc_calls WHERE list_contains(doms,'DNA_pol_A')
)
SELECT LOCUSTAG, gene, count(distinct regexp_replace(protein_id,'-T[0-9]+$','')) AS n FROM u GROUP BY 1,2;
SELECT gene, count(*) genomes, round(avg(n),2) FROM kx2_kc_gene_counts GROUP BY 1 ORDER BY 1;
