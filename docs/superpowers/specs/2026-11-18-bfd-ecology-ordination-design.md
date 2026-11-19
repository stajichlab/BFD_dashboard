# Design: ecology (FungalTraits) in the summary ordinations

Date: 2026-11-18. Status: draft for the user's approval. Approach A (a separate step) and the zoom-in rule were agreed in conversation; the rest below is my proposal.

## Intent
The summary site has PCoA ordinations of codon usage, amino-acid composition and gene structure, coloured by phylum and class. The user wants the same plots with ecology added as a colour variable, and zoomed-in ordinations for every clade that meets a threshold, because it is not known which clades are informative. The user knows taxonomy explains more than ecology. Within clades the question is still open.

## Facts (measured 2026-11-18 on the scratch build `db/BFD_traits_dev.duckdb`)
- Ordination unit: one genome per genus (`pick_complete_reps` in `bin/summary/summary_ordination.py`). FungalTraits is a genus label, so each point has exactly one label.
- 1,669 genus representatives (`is_genus_rep`); 1,492 (89.4%) have `primary_lifestyle`; 1,660 have `growth_form`. The ordination reps may differ slightly because they require complete feature data.
- 28 lifestyle values. Largest: plant_pathogen 259, wood_saprotroph 230, lichenized 185, soil_saprotroph 170, litter_saprotroph 121, animal_parasite 105, ectomycorrhizal 102, unspecified_saprotroph 87. Many have under 10 genera.
- Clades with at least 30 genus reps and at least 2 lifestyles with at least 5 genera: 3 phyla, 8 classes (Agaricomycetes 300, Sordariomycetes 236, Dothideomycetes 200, Leotiomycetes 84, Eurotiomycetes 66, Pichiomycetes 40, Tremellomycetes 36, Chytridiomycetes 35), 6 orders.
- `genome_all.parquet` has no ecology columns and is shared with the `bfd-build` session. `make_figures.py` `main()` and `figures_release.json` are also `bfd-build`'s.

## Goals
1. Global ordinations (same genus reps, same distances, same PCoA as `run_ordination`) coloured by `primary_lifestyle` and by `growth_form`, for each of codon, aa and structure.
2. A PERMANOVA of lifestyle on the whole tree, unrestricted and with lifestyle labels permuted only within order.
3. For every clade at phylum, class and order level that passes the threshold: its own PCoA (recomputed on that clade's genera), a scatter coloured by lifestyle, and a lifestyle PERMANOVA, for each feature set. Benjamini-Hochberg correction across clades within each level and feature set.
4. One new Quarto page, `17_ecology.qmd`, with the global figures, the test tables and the clade table.
5. Counts always reconcile: genera coloured + genera without a label = genera in the ordination.

## Non-goals and anti-goals
- No edit to `make_figures.py`, `_quarto.yml`, `06_data_notes.qmd`, `figure_index.py`, `figures_release.json`, `genome_all.parquet`, or `tf_*` files. No release upload from this work.
- Not in this design: distributions of genome traits by growth form (dimorphic vs yeast vs filamentous). That is a later, separate design.
- Not a claim that ecology drives genome content. Labels are genus-level and coarse, and lifestyle is confounded with ancestry. A significant result inside a clade is a pattern, not a cause.
- A failure case: a clade figure whose points are mostly unlabeled and look like a result.

## Approaches
1. **Separate step `summary_ecology.py`, run after `make_figures`.** Chosen by the user. It recomputes the global PCoA with the same code and seed, so coordinates match the existing ordination tables. It needs no change to shared files.
2. Add ecology inside `run_ordination`. Rejected: needs ecology columns in the shared `genome_all.parquet`.
3. Read the FungalTraits CSV in the plotting code. Rejected: duplicates the join now in the DB.

## Design
- **Entry point:** `python bin/summary/summary_ecology.py --tables TABLES_DIR --db BFD.duckdb --out OUT_DIR [--min-genera 30 --min-group 5 --levels phylum class order]`. Runs from `workflows/summary_ecology.sh` (SLURM). Reads `TABLES_DIR/{genome_all,codon_genome,aa_genome,gene_structure_genome}.parquet` read-only and `v_genome_labels` from the DB.
- **Traits:** `load_traits(db) -> DataFrame[species_prefix, primary_lifestyle, growth_form]` from `v_genome_labels` (`LOCUSTAG AS species_prefix`). Genera with a NULL label are drawn grey as "no label" and are not part of any test.
- **Reuse:** import `pick_complete_reps`, the matrix builders and the PCoA/distance code from `summary_ordination.py` and `lib/ordination.py`; `Out` and `assign_groups` from `make_figures.py`. No change to those modules.
- **Colouring:** lifestyle groups with fewer than `--min-group` genera in the plotted set are pooled as "other (n groups)". Legend outside the axes, two columns, 6 pt, as in PR #19.
- **Within a clade:** the feature matrix and the structure z-scoring are recomputed on the clade's genera; the zoom is its own ordination, not a crop of the global plot.
- **Tests:** PERMANOVA pseudo-F, R2, p with 999 permutations and seed 1 (as in `summary_ordination.py`). Lifestyles with fewer than `--min-group` genera in the tested set are dropped from the test and counted. Stratified permutation: labels are permuted within order; orders with one genus contribute nothing. Existing `caveat` column convention.
- **Outputs** under `OUT_DIR`:
  - `figures/ecology/E01_<matrix>_pcoa_lifestyle`, `E02_<matrix>_pcoa_growthform`, `E03_<level>_<clade>_<matrix>_pcoa_lifestyle` (svg and png).
  - `tables/ecology/ecology_traits_coverage.tsv`, `ecology_permanova_global.tsv`, `ecology_permanova_clades.tsv` (level, clade, matrix, n_genera, n_labeled, n_groups, pseudo_F, R2, p, q_BH, caveat), `ecology_clades_skipped.tsv` (clade, reason).
- **Page:** `docs/quarto_summary/17_ecology.qmd`, using the existing `_helpers.py`. Navbar entry is added after `bfd-build`'s PR #24 merges.
- **Dependencies:** needs the trait columns in the DB that the step reads. Development uses `db/BFD_traits_dev.duckdb`; the final run waits for the held production build (job 32195607) or the DB the user chooses.

## Acceptance
- The global PCoA coordinates equal `tables/ordination/<matrix>_coords.tsv` to 1e-9.
- For every figure, coloured + unlabeled = n genera plotted.
- Tests on a synthetic table: pooling, stratified permutation (labels never cross an order), BH correction, a clade below threshold goes to the skipped table with a reason.
- About 17 clades times 3 feature sets at the default threshold; the real count is printed.

## Left to the builder
- How `17_ecology.qmd` lays out ~51 clade figures (a table with links, or collapsible sections).
- Colour palette for up to 9 lifestyle groups plus "other" and "no label"; it must differ from the phylum palette.
- Whether `figure_index.py` picks up `figures/ecology` automatically; if it does not, ask `bfd-build` and do not edit it.

## Later (not part of this build)
- Dimorphic vs yeast vs filamentous: genome-trait distributions by `growth_form`, as its own design.
- Phylogenetic eigenvectors to separate ecology from ancestry (noted in `docs/FUNCTIONAL_ANALYSES.md`).
