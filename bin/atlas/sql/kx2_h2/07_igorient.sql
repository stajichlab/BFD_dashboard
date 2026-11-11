-- 07_igorient.sql
-- Purpose: per-genome, per-orientation intergenic distance statistics, then genus medians -> lineage medians.
-- Inputs: bfd_scratch.kx2_h2_ig, kx2_h2_genomes. Output: bfd_scratch.kx2_h2_igorient (+ a result query).
-- Filters: result query uses Q genomes, plus all Harpellales/Dimargaritales genomes (no Q genomes exist for them).
-- Mode = most frequent 25-bp bin start. f_ov = distance <= 0 (0 everywhere).
-- Report numbers: "ORIENTATION-SPECIFIC INTERGENIC" table.
CREATE OR REPLACE TABLE kx2_h2_igorient AS
SELECT i.genome, g.lineage, g.GENUS, g.q, g.busco_complete, g.N50_bp, i.orient, count(*) n, median(distance) med, mode(distance//25*25) mode25, avg((distance<100)::int) f_lt100, avg((distance<=0)::int) f_ov, avg((distance<=20)::int) f_le20
FROM kx2_h2_ig i JOIN kx2_h2_genomes g ON g.LOCUSTAG=i.genome WHERE g.lineage IS NOT NULL GROUP BY ALL;

-- result query
-- WITH gm AS (SELECT lineage, GENUS, orient, median(med) med, median(mode25) mode25, median(f_lt100) f_lt100, median(f_ov) f_ov, median(n) n FROM kx2_h2_igorient WHERE q OR lineage IN ('Harpellales','Dimargaritales') GROUP BY ALL)
-- SELECT lineage, orient, count(*) ng, round(median(med)) med, round(median(mode25)) mode25, round(median(f_lt100),3) lt100, round(median(f_ov),4) ov FROM gm GROUP BY ALL ORDER BY lineage, orient;
