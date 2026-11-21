-- REF_DIR=/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/config/refs
-- dashboard_views.sql: the dashboard's analysis layer, written into its own small file db/BFD_labels.duckdb.
-- Built by workflows/build_labels_db.sh. Run with the labels file as the current database and the pipeline
-- database attached read-only as `bfd`:
--     ATTACH '<pipeline BFD.duckdb>' AS bfd (READ_ONLY);   -- then run this file
-- Line 1 is replaced with the real config/refs path; @SOURCE_DB@ and @SOURCE_MTIME@ in labels_meta are replaced
-- by the build script. Idempotent.
--
-- Why a separate file: the pipeline owns build FACTS (views over parquet). This file holds analysis JUDGEMENTS
-- (quality thresholds, representatives, FungalTraits labels, interim SUBPHYLUM patches). They change often and must
-- not need a pipeline release. Results are stored as tables (about 24,000 rows), not views, because a DuckDB view
-- cannot reliably point at another attached file and the labels are cheap to store.
-- v_genome_traits and v_genome_labels keep their old names and columns, so queries written for them still work.
-- The rules are those of sql/bfd_traits.sql and sql/bfd_extra_views.sql (v_genome_labels), unchanged.

CREATE OR REPLACE TABLE ref_fungaltraits AS
SELECT GENUS,
       Growth_form_template AS growth_form_ft,
       primary_lifestyle,
       Secondary_lifestyle AS secondary_lifestyle,
       Plant_pathogenic_capacity_template AS plant_pathogenic_capacity,
       Decay_substrate_template AS decay_substrate,
       Animal_biotrophic_capacity_template AS animal_biotrophic_capacity
FROM (SELECT *, row_number() OVER (PARTITION BY GENUS ORDER BY jrk_template) AS rn
      FROM read_csv('/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/config/refs/fungaltraits_1.2.utf8.csv', header=true, all_varchar=true))
WHERE rn = 1 AND GENUS IS NOT NULL;
COMMENT ON TABLE ref_fungaltraits IS 'FungalTraits 1.2, one row per GENUS (first row by jrk_template when a genus is listed twice). Source: config/refs/fungaltraits_1.2.utf8.csv.';

CREATE OR REPLACE TABLE ref_growth_form_override AS
SELECT GENUS, max(nullif(growth_form, '')) AS growth_form_ovr
FROM read_csv('/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/config/refs/growth_form_overrides.tsv', header=true, delim='\t', all_varchar=true)
WHERE GENUS IS NOT NULL
GROUP BY GENUS;
COMMENT ON TABLE ref_growth_form_override IS 'Growth form for genera that FungalTraits lacks (config/refs/growth_form_overrides.tsv). An empty growth_form is ignored.';

CREATE OR REPLACE TABLE ref_class_subphylum AS
SELECT CLASS, min(SUBPHYLUM) AS SUBPHYLUM, min(source) AS source, min(note) AS note
FROM read_csv('/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/config/refs/class_subphylum.tsv', header=true, delim='\t', all_varchar=true)
WHERE CLASS IS NOT NULL
GROUP BY CLASS;
COMMENT ON TABLE ref_class_subphylum IS 'CLASS to SUBPHYLUM (config/refs/class_subphylum.tsv). source: taxonomy_py, atlas, or proposed (not yet checked). Classes not listed have no subphylum here.';

CREATE OR REPLACE TABLE ref_phylum_subphylum AS
SELECT PHYLUM, min(SUBPHYLUM) AS SUBPHYLUM, min(source) AS source, min(note) AS note
FROM read_csv('/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/config/refs/phylum_subphylum.tsv', header=true, delim='\t', all_varchar=true)
WHERE PHYLUM IS NOT NULL
GROUP BY PHYLUM;
COMMENT ON TABLE ref_phylum_subphylum IS 'PHYLUM to SUBPHYLUM fallback for genomes whose CLASS is missing or not in ref_class_subphylum (config/refs/phylum_subphylum.tsv).';

CREATE OR REPLACE TABLE ref_species_subphylum AS
SELECT SPECIES, min(SUBPHYLUM) AS SUBPHYLUM, min(source) AS source, min(note) AS note
FROM read_csv('/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/config/refs/species_subphylum.tsv', header=true, delim='\t', all_varchar=true)
WHERE SPECIES IS NOT NULL
GROUP BY SPECIES;
COMMENT ON TABLE ref_species_subphylum IS 'SPECIES name to SUBPHYLUM for genomes with no usable CLASS (config/refs/species_subphylum.tsv). Placeholder names until the taxonomy is fixed upstream. Beats the class and phylum tables, not an NCBI value.';

CREATE OR REPLACE TABLE v_genome_traits AS
SELECT t.LOCUSTAG, t.growth_form, t.growth_form_source, t.primary_lifestyle, t.secondary_lifestyle,
       t.plant_pathogenic_capacity, t.decay_substrate, t.animal_biotrophic_capacity, t.in_fungaltraits,
       CASE WHEN t.growth_form = 'yeast' AND t.CLASS IN ('Saccharomycetes','Pichiomycetes','Dipodascomycetes','Lipomycetes','Trigonopsidomycetes','Alloascoideomycetes','Sporopachydermiomycetes') THEN 'Saccharomycotina'
            WHEN t.growth_form = 'yeast' AND t.CLASS IN ('Schizosaccharomycetes','Taphrinomycetes','Pneumocystomycetes','Neolectomycetes') THEN 'Taphrinomycotina'
            WHEN t.growth_form = 'yeast' AND t.CLASS = 'Tremellomycetes' THEN 'Tremellomycetes'
            WHEN t.growth_form = 'yeast' AND t.CLASS IN ('Microbotryomycetes','Cystobasidiomycetes','Agaricostilbomycetes','Mixiomycetes','Spiculogloeomycetes','Classiculomycetes','Tritirachiomycetes') THEN 'Pucciniomycotina-yeast'
            WHEN t.growth_form = 'yeast' AND t.CLASS IN ('Malasseziomycetes','Ustilaginomycetes','Exobasidiomycetes','Moniliellomycetes') THEN 'Ustilaginomycotina-yeast'
            WHEN t.growth_form = 'yeast' THEN 'yeast-other'
            WHEN t.growth_form = 'dimorphic_yeast' AND t.CLASS IN ('Eurotiomycetes','Dothideomycetes','Sordariomycetes','Leotiomycetes','Lecanoromycetes','Pezizomycetes','Orbiliomycetes','Arthoniomycetes','Xylonomycetes','Geoglossomycetes','Lichinomycetes','Laboulbeniomycetes','Coniocybomycetes','Candelariomycetes') THEN 'Pezizomycotina-dimorphic'
            WHEN t.growth_form = 'dimorphic_yeast' AND t.CLASS IN ('Malasseziomycetes','Ustilaginomycetes','Exobasidiomycetes','Moniliellomycetes') THEN 'Ustilaginomycotina-dimorphic'
            WHEN t.growth_form = 'dimorphic_yeast' THEN 'dimorphic-other'
            WHEN t.growth_form = 'filamentous_mycelium' THEN 'filamentous'
            ELSE 'other' END AS yeast_origin,
       t.SUBPHYLUM, t.subphylum_source
FROM (
  SELECT s.LOCUSTAG, s.CLASS,
         coalesce(o.growth_form_ovr, f.growth_form_ft) AS growth_form,
         CASE WHEN o.growth_form_ovr IS NOT NULL THEN 'override'
              WHEN f.growth_form_ft IS NOT NULL THEN 'fungaltraits' END AS growth_form_source,
         f.primary_lifestyle, f.secondary_lifestyle, f.plant_pathogenic_capacity, f.decay_substrate,
         f.animal_biotrophic_capacity,
         (f.GENUS IS NOT NULL) AS in_fungaltraits,
         coalesce(s.SUBPHYLUM, sp.SUBPHYLUM, c.SUBPHYLUM, ph.SUBPHYLUM) AS SUBPHYLUM,
         CASE WHEN s.SUBPHYLUM IS NOT NULL THEN 'species'
              WHEN sp.SUBPHYLUM IS NOT NULL THEN 'species_table'
              WHEN c.SUBPHYLUM IS NOT NULL THEN 'class_table'
              WHEN ph.SUBPHYLUM IS NOT NULL THEN 'phylum_table' END AS subphylum_source
  FROM bfd.species s
  LEFT JOIN ref_fungaltraits f ON s.GENUS = f.GENUS
  LEFT JOIN ref_growth_form_override o ON s.GENUS = o.GENUS
  LEFT JOIN ref_species_subphylum sp ON s.SPECIES = sp.SPECIES
  LEFT JOIN ref_class_subphylum c ON s.CLASS = c.CLASS
  LEFT JOIN ref_phylum_subphylum ph ON s.PHYLUM = ph.PHYLUM
) t;
COMMENT ON TABLE v_genome_traits IS 'One row per genome (key LOCUSTAG): FungalTraits growth form and lifestyle by GENUS (in_fungaltraits false = genus not in FungalTraits 1.2), yeast_origin, and SUBPHYLUM (species.SUBPHYLUM if filled, else the species-name table, else the class table, else the phylum table; NULL when none). Same names and rules as atlas genome_labels.parquet.';

CREATE OR REPLACE TABLE v_genome_labels AS
WITH base AS (
  SELECT s.LOCUSTAG, s.ASMID, s.PHYLUM, s.CLASS, s."ORDER" AS ORDER_, s.FAMILY, s.GENUS, s.SPECIES, s.STRAIN,
         s.NCBI_TAXONID, sa.TRANSL_TABLE AS transl_table,
         a.annotation_status, a.TOTAL_LENGTH AS total_length_bp, a.contig_count, a.N50_bp,
         a.GC_PERCENT AS gc_pct, a.masked_pct,
         b.complete_pct AS busco_complete, b.single_pct AS busco_single, b.duplicated_pct AS busco_duplicated,
         b.fragmented_pct AS busco_fragmented, b.missing_pct AS busco_missing,
         t.* EXCLUDE (LOCUSTAG)
  FROM bfd.species s
  LEFT JOIN v_genome_traits t ON s.LOCUSTAG = t.LOCUSTAG
  LEFT JOIN bfd.samples sa ON s.LOCUSTAG = sa.LOCUSTAG
  LEFT JOIN bfd.asm_stats a ON s.ASMID = a.ASMID
  LEFT JOIN bfd.busco_genome b ON s.ASMID = b.ASMID
), cls_thr AS (
  SELECT CLASS, quantile_cont(busco_complete, 0.25) AS p25, quantile_cont(busco_complete, 0.5) AS p50
  FROM base WHERE annotation_status = 'annotated' AND CLASS IS NOT NULL
  GROUP BY CLASS HAVING count(*) >= 20
), phy_thr AS (
  SELECT PHYLUM, quantile_cont(busco_complete, 0.25) AS p25, quantile_cont(busco_complete, 0.5) AS p50
  FROM base WHERE annotation_status = 'annotated' AND PHYLUM IS NOT NULL
  GROUP BY PHYLUM
), flagged AS (
  SELECT l.*,
         coalesce(c.p25, p.p25) AS busco_thr_p25,
         coalesce(c.p50, p.p50) AS busco_thr_p50,
         CASE WHEN c.CLASS IS NOT NULL THEN 'class' WHEN p.PHYLUM IS NOT NULL THEN 'phylum' END AS busco_thr_level,
         (l.annotation_status = 'annotated') AS q_ann,
         (l.busco_complete >= coalesce(c.p25, p.p25)) AS q_busco_rel,
         (l.busco_complete >= coalesce(c.p50, p.p50)) AS q_busco_rel50,
         (l.busco_duplicated <= 10) AS q_dup,
         (l.N50_bp >= 50000) AS q_n50,
         (row_number() OVER (PARTITION BY l.SPECIES
                             ORDER BY l.busco_complete DESC NULLS LAST, l.N50_bp DESC NULLS LAST, l.LOCUSTAG) = 1
          AND l.annotation_status = 'annotated' AND l.SPECIES IS NOT NULL) AS is_representative
  FROM base l
  LEFT JOIN cls_thr c ON l.CLASS = c.CLASS
  LEFT JOIN phy_thr p ON l.PHYLUM = p.PHYLUM
)
SELECT f.*,
       (f.is_representative AND f.GENUS IS NOT NULL AND
        row_number() OVER (PARTITION BY f.GENUS, f.is_representative
                           ORDER BY f.busco_complete DESC NULLS LAST, f.N50_bp DESC NULLS LAST, f.LOCUSTAG) = 1) AS is_genus_rep
FROM flagged f;
COMMENT ON TABLE v_genome_labels IS 'One row per genome: taxonomy (ORDER_ has a trailing underscore), assembly and BUSCO values, quality flags q_ann, q_busco_rel, q_busco_rel50, q_dup, q_n50, and is_representative (one genome per SPECIES) and is_genus_rep (one per GENUS). Same rules and column names as atlas genome_labels.parquet (bin/atlas/build_labels.py); FungalTraits columns (growth_form, growth_form_source, primary_lifestyle, secondary_lifestyle, plant_pathogenic_capacity, decay_substrate, animal_biotrophic_capacity, in_fungaltraits, yeast_origin) are a GENUS-level label from FungalTraits 1.2 (in_fungaltraits false = genus not listed); SUBPHYLUM is NULL when the class has no entry in ref_class_subphylum (subphylum_source says where it came from). Use WHERE is_representative before counting genomes per clade, or strains are counted as species. Key LOCUSTAG.';

CREATE OR REPLACE TABLE labels_meta AS
SELECT '@SOURCE_DB@' AS source_db, CAST(@SOURCE_MTIME@ AS BIGINT) AS source_mtime_epoch, current_timestamp AS built_at,
       (SELECT count(*) FROM v_genome_labels) AS n_genomes;
COMMENT ON TABLE labels_meta IS 'Which pipeline database these labels were built from (path and file mtime) and when. bfd_connect warns when the pipeline file is newer.';
