-- BUSCO gap by group: genome-mode vs protein-mode BUSCO (fungi_odb12), as a check on the funannotate annotation.
-- Inputs: bfd.species, bfd.busco_genome (genome-mode complete %), bfd.asm_stats, and bfd_scratch.kx2_pbusco, which holds the
-- protein-mode values from protein_busco_169.tsv (columns locustag, p_c, p_s, p_d, p_f, p_m; load with
--   CREATE TABLE bfd_scratch.kx2_pbusco AS SELECT * FROM read_csv('protein_busco_169.tsv', delim='\t', header=true);
-- renaming the columns to locustag, p_c, p_s, p_d, p_f, p_m).
-- gap = genome-mode complete % minus protein-mode complete %. A large positive gap (the BFD team's "BUSCO gap" set used
-- genome >= 90% and protein < 50%) means the assembly contains genes that the annotation did not recover.
-- Output: one row per group (run 2026-11-11). Result for Coemansia (107 strains): genome 80.4, protein 82.4, median gap -0.2,
-- range -3.2 to +3.9, 0 strains with a gap above 10 points.
WITH j AS (
  SELECT s.LOCUSTAG, s.GENUS, s."ORDER" AS ord, s.PHYLUM, b.complete_pct AS g_c, p.p_c, a.N50_bp, (b.complete_pct - p.p_c) AS gap,
    CASE WHEN s.GENUS='Coemansia' THEN '1 Coemansia (strains)'
         WHEN s."ORDER"='Kickxellales' THEN '2 other Kickxellales'
         WHEN s."ORDER"='Harpellales' THEN '3 Harpellales'
         WHEN s."ORDER" IN ('Dimargaritales','Ramicandelaberales') THEN '4 Dimargaritales+Ramicandelaber'
         WHEN s."ORDER"='Zoopagales' THEN '5 Zoopagales'
         WHEN s."ORDER"='Entomophthorales' OR s.PHYLUM='Basidiobolomycota' THEN '6 Entomophthorales+Basidiobolus'
         WHEN s.PHYLUM IN ('Mucoromycota','Mortierellomycota','Glomeromycota') THEN '7 Mucoromycota group'
         ELSE '8 root (chytrids etc.)' END AS grp
  FROM bfd_scratch.kx2_pbusco p
  JOIN species s ON s.LOCUSTAG = p.locustag
  JOIN busco_genome b ON b.ASMID = s.ASMID
  JOIN asm_stats a ON a.ASMID = s.ASMID)
SELECT grp, count(*) AS n, round(median(g_c),1) AS genome_busco, round(median(p_c),1) AS protein_busco,
       round(median(gap),1) AS median_gap, round(min(gap),1) AS min_gap, round(max(gap),1) AS max_gap,
       count(*) FILTER (WHERE gap > 10) AS n_gap_gt10, count(*) FILTER (WHERE gap < -5) AS n_protein_higher
FROM j GROUP BY grp ORDER BY grp;
