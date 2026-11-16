-- 02_kx_genomes.sql  (J12 rebuild; RECONSTRUCTED, see below)
-- Purpose: species representatives of the Kickxellomycota study set and comparison lineages, with a hand-defined group label grp.
-- STATUS: the original definition was only described in prose (kickxellomycota/README.md: "species JOIN asm_stats JOIN busco_genome,
--   row_number() OVER (PARTITION BY SPECIES ORDER BY complete_pct DESC, N50_bp DESC, LOCUSTAG) = 1, CASE on PHYLUM/CLASS/ORDER for grp").
--   The CASE and the genus-level thinning of large groups below are reconstructed from the MotherDuck table kx_genomes (group names,
--   counts and per-group LOCUSTAG hashes). J12_rebuild.md reports how closely the rebuilt table matches it.
-- Thinning (reconstructed): groups where the MotherDuck table has exactly one genome per genus (six Basidiomycota classes, six
--   Pezizomycotina classes, Sacch: Dipodascomycetes, Pichiomycetes, Saccharomycetes) keep only the best genome per genus.
-- Output: bfd_scratch.kx_genomes (MotherDuck: 1,919 rows)
CREATE OR REPLACE TABLE bfd_scratch.kx_genomes AS
WITH g AS (
  SELECT s.ASMID, s.LOCUSTAG, s.PHYLUM, s.CLASS, s.ORDER_, s.FAMILY, s.GENUS, s.SPECIES, s.STRAIN,
         a.total_length_bp, a.N50_bp, a.contig_count, a.gc_pct, a.masked_pct, a.annotation_status,
         b.complete_pct AS busco_complete, b.duplicated_pct AS busco_dup, b.fragmented_pct AS busco_frag, b.missing_pct AS busco_missing
  FROM bfd_scratch.j12_species s
  JOIN bfd_scratch.j12_asm_stats a USING (ASMID)
  LEFT JOIN bfd.busco_genome b USING (ASMID)
  WHERE a.annotation_status = 'annotated'),
r AS (
  SELECT *, row_number() OVER (PARTITION BY SPECIES ORDER BY busco_complete DESC NULLS LAST, N50_bp DESC NULLS LAST, LOCUSTAG) AS sp_rk
  FROM g),
lab AS (
  SELECT *,
    CASE
      WHEN ORDER_ = 'Kickxellales' AND GENUS = 'Coemansia' THEN 'Kx:Coemansia'
      WHEN ORDER_ = 'Kickxellales' AND GENUS IN ('Kickxella','Linderina','Dipsacomyces','Martensiomyces') THEN 'Kx:Kickxellales-core-other'
      WHEN ORDER_ = 'Kickxellales' THEN 'Kx:Kickxellales-early'
      WHEN ORDER_ = 'Harpellales' THEN 'Kx:Harpellales'
      WHEN ORDER_ = 'Dimargaritales' THEN 'Kx:Dimargaritales'
      WHEN ORDER_ = 'Ramicandelaberales' THEN 'Kx:Ramicandelaber'
      WHEN ORDER_ = 'Zoopagales' THEN 'Zoopagales'
      WHEN ORDER_ = 'Mucorales' THEN 'Mucorales'
      WHEN ORDER_ = 'Endogonales' THEN 'Endogonales'
      WHEN ORDER_ = 'Umbelopsidales' THEN 'Umbelopsidales'
      WHEN PHYLUM = 'Mortierellomycota' THEN 'Mortierellales'
      WHEN PHYLUM = 'Glomeromycota' THEN 'Glomeromycota'
      WHEN PHYLUM = 'Basidiobolomycota' THEN 'Basidiobolus'
      WHEN PHYLUM = 'Entomophthoromycota' THEN 'Entomophthorales'
      WHEN PHYLUM = 'Blastocladiomycota' THEN 'Blastocladiomycota'
      WHEN PHYLUM = 'Neocallimastigomycota' THEN 'Neocallimastigomycota'
      WHEN PHYLUM IN ('Chytridiomycota','Monoblepharomycota') THEN 'Chytrid+Monobleph'
      WHEN PHYLUM IN ('Rozellomycota','Aphelidiomycota') THEN 'Rozello+Aphelid'
      WHEN PHYLUM = 'Microsporidia' THEN 'Microsporidia'
      WHEN CLASS IN ('Alloascoideomycetes','Dipodascomycetes','Lipomycetes','Pichiomycetes','Saccharomycetes','Sporopachydermiomycetes','Trigonopsidomycetes') THEN 'Sacch:' || CLASS
      WHEN CLASS = 'Schizosaccharomycetes' THEN 'Schizosaccharomyces'
      WHEN CLASS = 'Taphrinomycetes' THEN 'Taphrinomycetes'
      WHEN CLASS = 'Pneumocystomycetes' THEN 'Pneumocystis'
      WHEN CLASS IN ('Neolectomycetes','Saitoellomycetes','Novakomycetes') THEN 'Taphrino-other'
      WHEN CLASS IN ('Dothideomycetes','Eurotiomycetes','Leotiomycetes','Orbiliomycetes','Pezizomycetes','Sordariomycetes') THEN 'Pez:' || CLASS
      WHEN CLASS IN ('Agaricomycetes','Exobasidiomycetes','Microbotryomycetes','Pucciniomycetes','Tremellomycetes','Ustilaginomycetes') THEN CLASS
      WHEN CLASS = 'Malasseziomycetes' THEN 'Malassezia'
      WHEN CLASS = 'Wallemiomycetes' THEN 'Wallemia'
    END AS grp
  FROM r WHERE sp_rk = 1),
gen AS (
  SELECT *, row_number() OVER (PARTITION BY GENUS ORDER BY busco_complete DESC NULLS LAST, N50_bp DESC NULLS LAST, LOCUSTAG) AS gen_rk
  FROM lab WHERE grp IS NOT NULL)
SELECT ASMID, LOCUSTAG, PHYLUM, CLASS, ORDER_, FAMILY, GENUS, SPECIES, STRAIN, total_length_bp, N50_bp, contig_count, gc_pct, masked_pct,
       annotation_status, busco_complete, busco_dup, busco_frag, busco_missing, grp, TRUE AS is_sp_rep, (gen_rk = 1) AS is_genus_rep
FROM gen
WHERE gen_rk = 1
   OR grp NOT IN ('Agaricomycetes','Exobasidiomycetes','Microbotryomycetes','Pucciniomycetes','Tremellomycetes','Ustilaginomycetes',
                  'Pez:Dothideomycetes','Pez:Eurotiomycetes','Pez:Leotiomycetes','Pez:Orbiliomycetes','Pez:Pezizomycetes','Pez:Sordariomycetes',
                  'Sacch:Dipodascomycetes','Sacch:Pichiomycetes','Sacch:Saccharomycetes');
-- check: per-group count and md5 of the sorted LOCUSTAG list (compare with the MotherDuck table, see J12_rebuild.md)
SELECT grp, count(*) AS n, md5(string_agg(LOCUSTAG, ',' ORDER BY LOCUSTAG)) AS h FROM bfd_scratch.kx_genomes GROUP BY 1 ORDER BY 1;
