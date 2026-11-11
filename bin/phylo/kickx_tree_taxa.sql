-- Taxon selection for the Kickxellomycota phylogenies (backbone tree and Coemansia tree).
-- Run with bin/atlas/run_duckdb_job.py, pointing --atlas at a purgeable bench folder, e.g.
--   ~/projects/BFD_dashboard/bench/kickx_tree/   (outputs land in <atlas>/taxa/)
-- Placeholders: {T} = BFD parquet tables dir, {OUT} = block output path (see run_duckdb_job.py).
--
-- Representative rule (same as bfd_scratch.kw_rep_genomes): annotated genomes only; within a genus keep the
-- genome with the highest BUSCO complete_pct, then the largest N50, then the lowest ASMID.
-- BUSCO complete_pct is fungi_odb12 here (busco_genome.lineage); it is used only to PICK genomes, not as a tree marker set.
--
-- Backbone design (the question: where did genome size, intron density and gene content change?):
--   ingroup     all Kickxellomycota genera, no BUSCO floor (several are single-genome and fragmented)
--   sister      Zoopagomycota, Entomophthoromycota, Basidiobolomycota: all genera with BUSCO >= 60
--   mucoromycot up to 8 genera of Mucorales, up to 2 genera per other order of Mucoromycota, Mortierellomycota,
--               Glomeromycota, with BUSCO >= 85 (Mucorales is the baseline used in the Analysis E accounting)
--   root        one genus per order of Blastocladiomycota, Chytridiomycota, Monoblepharomycota,
--               Neocallimastigomycota with BUSCO >= 75
-- The zygomycete backbone order is contested, so the tree must be run and compared on at least two backbones
-- (concatenation with a site-heterogeneous model, and a gene-tree summary method); see FUNCTIONAL_ANALYSES.md, Analysis E.
--
-- Coemansia design (the question: relationships among strains, and how many independent lineages there are):
--   all annotated Coemansia genomes with BUSCO >= 60, plus the other core Kickxellales genera, Dimargaritales and
--   Harpellales as outgroups. Strains are NOT species: 129 of 146 Coemansia representatives are unnamed "sp.".
--   n50_ge_50kb marks the high-quality assemblies (BUSCO 87-91%) that anchor the tree.

CREATE OR REPLACE VIEW g AS
SELECT s.ASMID, s.LOCUSTAG, s.PHYLUM, s."CLASS" AS cls, s."ORDER" AS ord, s.FAMILY, s.GENUS, s.SPECIES, s.STRAIN,
       b.complete_pct, a.N50_bp, a.total_length_bp, a.annotation_status
FROM '{T}/species.parquet' s
JOIN '{T}/busco_genome.parquet' b USING (ASMID)
JOIN '{T}/asm_stats.parquet' a USING (ASMID)
WHERE s.LOCUSTAG IS NOT NULL AND a.annotation_status = 'annotated';

CREATE OR REPLACE VIEW genus_rep AS
SELECT * FROM (
  SELECT *, row_number() OVER (PARTITION BY GENUS ORDER BY complete_pct DESC NULLS LAST, N50_bp DESC NULLS LAST, ASMID) AS rn
  FROM g WHERE GENUS IS NOT NULL
) WHERE rn = 1;

CREATE OR REPLACE VIEW genus_rep_ranked AS
SELECT *, row_number() OVER (PARTITION BY PHYLUM, ord ORDER BY complete_pct DESC NULLS LAST, N50_bp DESC NULLS LAST, ASMID) AS ord_rank
FROM genus_rep;

-- @output taxa/backbone_taxa.tsv
COPY (
  SELECT 'ingroup' AS role, LOCUSTAG, ASMID, PHYLUM, cls, ord, FAMILY, GENUS, SPECIES, complete_pct, N50_bp
  FROM genus_rep_ranked WHERE PHYLUM = 'Kickxellomycota'
  UNION ALL
  SELECT 'sister', LOCUSTAG, ASMID, PHYLUM, cls, ord, FAMILY, GENUS, SPECIES, complete_pct, N50_bp
  FROM genus_rep_ranked WHERE PHYLUM IN ('Zoopagomycota','Entomophthoromycota','Basidiobolomycota') AND complete_pct >= 60
  UNION ALL
  SELECT 'mucoromycot', LOCUSTAG, ASMID, PHYLUM, cls, ord, FAMILY, GENUS, SPECIES, complete_pct, N50_bp
  FROM genus_rep_ranked
  WHERE PHYLUM IN ('Mucoromycota','Mortierellomycota','Glomeromycota') AND complete_pct >= 85
    AND ord_rank <= CASE WHEN ord = 'Mucorales' THEN 8 ELSE 2 END
  UNION ALL
  SELECT 'root', LOCUSTAG, ASMID, PHYLUM, cls, ord, FAMILY, GENUS, SPECIES, complete_pct, N50_bp
  FROM genus_rep_ranked
  WHERE PHYLUM IN ('Blastocladiomycota','Chytridiomycota','Monoblepharomycota','Neocallimastigomycota')
    AND complete_pct >= 75 AND ord_rank = 1
  ORDER BY 1, 4, 6, 8
) TO '{OUT}' (FORMAT CSV, DELIMITER '\t', HEADER);

-- @output taxa/coemansia_taxa.tsv
COPY (
  SELECT 'ingroup' AS role, LOCUSTAG, ASMID, PHYLUM, cls, ord, FAMILY, GENUS, SPECIES, STRAIN, complete_pct, N50_bp,
         (N50_bp >= 50000) AS n50_ge_50kb
  FROM g WHERE GENUS = 'Coemansia' AND complete_pct >= 60
  UNION ALL
  SELECT 'outgroup', LOCUSTAG, ASMID, PHYLUM, cls, ord, FAMILY, GENUS, SPECIES, STRAIN, complete_pct, N50_bp,
         (N50_bp >= 50000)
  FROM genus_rep_ranked
  WHERE PHYLUM = 'Kickxellomycota' AND GENUS <> 'Coemansia'
  ORDER BY 1, 8, 11 DESC
) TO '{OUT}' (FORMAT CSV, DELIMITER '\t', HEADER);
