-- 01_genomes.sql
-- Purpose: H2 genome set: species representatives with hand-defined lineage labels and quality flags.
-- Inputs: bfd_scratch.kx_genomes (1,919 species reps with grp), bfd.asm_stats, joined on ASMID.
-- Output: bfd_scratch.kx2_h2_genomes
-- Definitions:
--   lineage: hand-defined mapping from kx_genomes.grp (CASE below). KxCore = Coemansia + Kickxellales-core-other
--     (Coemansia, Kickxella, Linderina, Dipsacomyces, Martensiomyces). Saccharomycotina = all grp 'Sacch:%'.
--     Pezizomycotina = all grp 'Pez:%'. Not assigned (NULL): other lineages, excluded downstream.
--   q  = busco_complete >= 85 AND N50_bp >= 50000 (Q filter). Q70 (busco>=70, N50>=25000) is computed in step 05.
-- Report numbers: group membership/genome counts for every table in the report.
-- Run against database bfd_scratch (write). Read-only on bfd.
CREATE OR REPLACE TABLE kx2_h2_genomes AS
SELECT g.LOCUSTAG, g.ASMID, g.GENUS, g.SPECIES, g.grp,
 CASE WHEN g.grp IN ('Kx:Coemansia','Kx:Kickxellales-core-other') THEN 'KxCore'
      WHEN g.grp='Kx:Harpellales' THEN 'Harpellales'
      WHEN g.grp='Kx:Dimargaritales' THEN 'Dimargaritales'
      WHEN g.grp='Zoopagales' THEN 'Zoopagales'
      WHEN g.grp='Mucorales' THEN 'Mucorales'
      WHEN g.grp LIKE 'Sacch:%' THEN 'Saccharomycotina'
      WHEN g.grp='Schizosaccharomyces' THEN 'Schizosaccharomyces'
      WHEN g.grp='Taphrinomycetes' THEN 'Taphrinomycetes'
      WHEN g.grp='Pneumocystis' THEN 'Pneumocystis'
      WHEN g.grp='Microsporidia' THEN 'Microsporidia'
      WHEN g.grp='Malassezia' THEN 'Malassezia'
      WHEN g.grp='Wallemia' THEN 'Wallemia'
      WHEN g.grp LIKE 'Pez:%' THEN 'Pezizomycotina'
      WHEN g.grp IN ('Kx:Kickxellales-early','Kx:Ramicandelaber') THEN 'KxOther'
 END AS lineage,
 g.busco_complete, g.N50_bp, g.contig_count, a.total_length_bp, a.masked_bases, a.total_n_bases, a.n_gap_count, a.t2t_scaffolds, a.telomere_fwd, a.telomere_rev,
 (g.busco_complete>=85 AND g.N50_bp>=50000) AS q
FROM kx_genomes g JOIN bfd.asm_stats a USING(ASMID);

-- check query run afterwards (counts per lineage, all vs Q)
-- SELECT lineage, count(*) n, count(distinct GENUS) g, sum(q::int) nq, count(distinct CASE WHEN q THEN GENUS END) gq FROM kx2_h2_genomes GROUP BY 1 ORDER BY 1;
