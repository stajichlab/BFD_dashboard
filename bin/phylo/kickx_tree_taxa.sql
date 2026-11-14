-- Taxon selection for the Kickxellomycota phylogenies (backbone tree and Coemansia tree).
-- Run with bin/atlas/run_duckdb_job.py, pointing --atlas at a purgeable bench folder, e.g.
--   ~/projects/BFD_dashboard/bench/kickx_tree/   (outputs land in <atlas>/taxa/)
-- Placeholders: {T} = BFD parquet tables dir, {OUT} = block output path (see run_duckdb_job.py).
--
-- Genome choice within a genus (revised 2026-11-11 after the first backbone run, see workflows/phylo/README.md):
--   rank by (annotation present: >= 3,000 proteins), then BUSCO SINGLE-copy complete (single_pct), then N50, then ASMID.
--   The first draft ranked by BUSCO complete_pct, which counts duplicated markers and so favoured duplicated
--   assemblies. In the first backbone run the genomes with BUSCO duplicated 52-70% (Rhizopus, Batkoa, Entomophthora,
--   Conidiobolus, Tieghemiomyces, Dispira, Polychytrium) recovered 14-426 of 1,089 fungi_odb12 markers in phyling align
--   against 600-1,000 for low-duplication genomes, and Massospora (BUSCO 71.7%) had only 415 proteins in BFD and 0 markers.
--   BUSCO complete_pct is fungi_odb12 (busco_genome.lineage); BUSCO is only used to PICK genomes, not as a tree marker set.
--
-- Backbone design (the question: where did genome size, intron density and gene content change?):
--   ingroup     all Kickxellomycota genera, no floor; flags low_busco (complete < 60), high_dup (duplicated > 10, the atlas
--               Q-dup threshold), low_prot (< 3,000 proteins) mark genomes whose marker recovery must be checked
--   sister      Zoopagomycota, Entomophthoromycota, Basidiobolomycota: all genera with >= 3,000 proteins and single_pct >= 40
--   mucoromycot up to 8 genera of Mucorales, up to 2 genera per other order of Mucoromycota, Mortierellomycota,
--               Glomeromycota, with >= 3,000 proteins and single_pct >= 60 (Mucorales is the Analysis E baseline)
--   root        one genus per order of Blastocladiomycota, Chytridiomycota, Monoblepharomycota,
--               Neocallimastigomycota with >= 3,000 proteins and single_pct >= 60
-- The zygomycete backbone order is contested, so the tree must be run and compared on at least two backbones
-- (concatenation with a site-heterogeneous model, and a gene-tree summary method); see FUNCTIONAL_ANALYSES.md, Analysis E.
--
-- Coemansia design (the question: relationships among strains, and how many independent lineages there are):
--   all annotated Coemansia genomes with BUSCO complete >= 60 and >= 3,000 proteins (all strains are kept, none are
--   dropped for duplication, but high_dup marks them), plus the other Kickxellomycota genera as outgroups chosen by
--   the same single-copy rule. Strains are NOT species: 129 of 146 Coemansia representatives are unnamed "sp.".
--   n50_ge_50kb marks the high-quality assemblies (BUSCO 87-91%) that anchor the tree.

CREATE OR REPLACE VIEW nprot AS
SELECT split_part(gene_id, '_', 1) AS LOCUSTAG, count(*) AS n_proteins
FROM '{T}/gene_proteins.parquet' GROUP BY 1;

CREATE OR REPLACE VIEW g AS
SELECT s.ASMID, s.LOCUSTAG, s.PHYLUM, s."CLASS" AS cls, s."ORDER" AS ord, s.FAMILY, s.GENUS, s.SPECIES, s.STRAIN,
       b.complete_pct, b.duplicated_pct, b.single_pct, a.N50_bp, a.total_length_bp, a.annotation_status,
       coalesce(n.n_proteins, 0) AS n_proteins
FROM '{T}/species.parquet' s
JOIN '{T}/busco_genome.parquet' b USING (ASMID)
JOIN '{T}/asm_stats.parquet' a USING (ASMID)
LEFT JOIN nprot n ON n.LOCUSTAG = s.LOCUSTAG
WHERE s.LOCUSTAG IS NOT NULL AND a.annotation_status = 'annotated';

CREATE OR REPLACE VIEW genus_rep AS
SELECT * FROM (
  SELECT *, row_number() OVER (
           PARTITION BY GENUS
           ORDER BY (n_proteins >= 3000) DESC, single_pct DESC NULLS LAST, N50_bp DESC NULLS LAST, ASMID) AS rn
  FROM g WHERE GENUS IS NOT NULL
) WHERE rn = 1;

-- order rank only among genomes that have a usable annotation, so a genus with an empty annotation never takes a slot
CREATE OR REPLACE VIEW genus_rep_ranked AS
SELECT *, row_number() OVER (PARTITION BY PHYLUM, ord ORDER BY single_pct DESC NULLS LAST, N50_bp DESC NULLS LAST, ASMID) AS ord_rank
FROM genus_rep WHERE n_proteins >= 3000;

-- @output taxa/backbone_taxa.tsv
COPY (
  SELECT 'ingroup' AS role, LOCUSTAG, ASMID, PHYLUM, cls, ord, FAMILY, GENUS, SPECIES, complete_pct, duplicated_pct, single_pct,
         n_proteins, N50_bp, (complete_pct < 60) AS low_busco, (duplicated_pct > 10) AS high_dup, (n_proteins < 3000) AS low_prot
  FROM genus_rep WHERE PHYLUM = 'Kickxellomycota'
  UNION ALL
  SELECT 'sister', LOCUSTAG, ASMID, PHYLUM, cls, ord, FAMILY, GENUS, SPECIES, complete_pct, duplicated_pct, single_pct,
         n_proteins, N50_bp, (complete_pct < 60), (duplicated_pct > 10), (n_proteins < 3000)
  FROM genus_rep_ranked WHERE PHYLUM IN ('Zoopagomycota','Entomophthoromycota','Basidiobolomycota') AND single_pct >= 40
  UNION ALL
  SELECT 'mucoromycot', LOCUSTAG, ASMID, PHYLUM, cls, ord, FAMILY, GENUS, SPECIES, complete_pct, duplicated_pct, single_pct,
         n_proteins, N50_bp, (complete_pct < 60), (duplicated_pct > 10), (n_proteins < 3000)
  FROM genus_rep_ranked
  WHERE PHYLUM IN ('Mucoromycota','Mortierellomycota','Glomeromycota') AND single_pct >= 60
    AND ord_rank <= CASE WHEN ord = 'Mucorales' THEN 8 ELSE 2 END
  UNION ALL
  SELECT 'root', LOCUSTAG, ASMID, PHYLUM, cls, ord, FAMILY, GENUS, SPECIES, complete_pct, duplicated_pct, single_pct,
         n_proteins, N50_bp, (complete_pct < 60), (duplicated_pct > 10), (n_proteins < 3000)
  FROM genus_rep_ranked
  WHERE PHYLUM IN ('Blastocladiomycota','Chytridiomycota','Monoblepharomycota','Neocallimastigomycota')
    AND single_pct >= 60 AND ord_rank = 1
  ORDER BY 1, 4, 6, 8
) TO '{OUT}' (FORMAT CSV, DELIMITER '\t', HEADER);

-- @output taxa/coemansia_taxa.tsv
COPY (
  SELECT 'ingroup' AS role, LOCUSTAG, ASMID, PHYLUM, cls, ord, FAMILY, GENUS, SPECIES, STRAIN, complete_pct, duplicated_pct,
         single_pct, n_proteins, N50_bp, (N50_bp >= 50000) AS n50_ge_50kb, (duplicated_pct > 10) AS high_dup
  FROM g WHERE GENUS = 'Coemansia' AND complete_pct >= 60 AND n_proteins >= 3000
  UNION ALL
  SELECT 'outgroup', LOCUSTAG, ASMID, PHYLUM, cls, ord, FAMILY, GENUS, SPECIES, STRAIN, complete_pct, duplicated_pct,
         single_pct, n_proteins, N50_bp, (N50_bp >= 50000), (duplicated_pct > 10)
  FROM genus_rep
  WHERE PHYLUM = 'Kickxellomycota' AND GENUS <> 'Coemansia'
  ORDER BY 1, 8, 11 DESC
) TO '{OUT}' (FORMAT CSV, DELIMITER '\t', HEADER);
