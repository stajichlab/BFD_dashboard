-- 00_kw_rep_genomes.sql  (J12 rebuild of an input table that the H1-H4 SQL reads)
-- Purpose: one representative annotated genome per species (best BUSCO complete, then N50, then ASMID).
-- Origin: kingdom_domains/README.md (the SQL shown there is not exactly what built the MotherDuck table; the MotherDuck table has the
--         12 columns below and no rank column, so this file writes those 12 columns).
-- Reads  : bfd.busco_genome and the J12 shim views bfd_scratch.j12_species / j12_asm_stats, which give the MotherDuck column names
--          (ORDER_, total_length_bp, gc_pct) that the local views name "ORDER", TOTAL_LENGTH, GC_PERCENT.
-- Output : bfd_scratch.kw_rep_genomes  (MotherDuck: 7,435 rows, 1,669 genera)
CREATE OR REPLACE TABLE bfd_scratch.kw_rep_genomes AS
WITH g AS (
  SELECT s.ASMID, s.LOCUSTAG, s.PHYLUM, s.CLASS, s.ORDER_, s.FAMILY, s.GENUS, s.SPECIES,
         b.complete_pct, b.duplicated_pct, a.N50_bp, a.total_length_bp,
         row_number() OVER (PARTITION BY s.SPECIES ORDER BY b.complete_pct DESC NULLS LAST, a.N50_bp DESC NULLS LAST, s.ASMID) AS rk
  FROM bfd_scratch.j12_species s
  LEFT JOIN bfd.busco_genome b USING (ASMID)
  LEFT JOIN bfd_scratch.j12_asm_stats a USING (ASMID)
  WHERE a.annotation_status = 'annotated')
SELECT * EXCLUDE (rk) FROM g WHERE rk = 1;
-- check
SELECT count(*) AS n_species_reps, count(DISTINCT GENUS) AS n_genera FROM bfd_scratch.kw_rep_genomes;
