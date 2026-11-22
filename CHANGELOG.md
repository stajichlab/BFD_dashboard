# Changelog

All notable changes to this project are documented in this file.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

The version number is bumped when a feature or fix is added. Until the project
reaches a stable 1.0, minor version bumps mean "something was added or changed".

## [Unreleased]

### Changed
- Switch to the pipeline database (2026-11-21). `lib/bfd_connect.py` now defaults to `Fungi_BFD_runs/BFD_function/db/BFD.duckdb` and, when that default is used and `db/BFD_labels.duckdb` exists, attaches the labels file automatically (`BFD_LABELS_DB=""` turns it off). `workflows/build_labels_db.sh` defaults to the pipeline copy. Checked against the rebuilt pipeline database (job 32361527, 57 views, species 24,363): the labels file (24,363 rows) equals the old dashboard `v_genome_labels` in every column except SUBPHYLUM for 173 genomes: 165 Microsporidia to Microsporidiomycotina (PR #32) and 8 Antarctolithica to Antarctolithicomycotina (PR #31). `workflows/build_bfd_db.sh` and `db/BFD.duckdb` are legacy and are not deleted.
- `bin/summary/summary_ecology.py` and `workflows/summary_ecology.sh` take `--labels` / `LABELS` (default `db/BFD_labels.duckdb`) and read `v_genome_labels` from it; the database default is the pipeline copy, which has no `v_genome_labels`.

### Added
- Dashboard labels layer, step 1 of moving the full DB build to the pipeline: `sql/dashboard_views.sql` and `workflows/build_labels_db.sh` build `db/BFD_labels.duckdb` (v_genome_labels, v_genome_traits and the ref_* tables as stored tables, plus `labels_meta` with the source database and its mtime) from an attached pipeline database. `lib/bfd_connect.py` takes `labels_db` or `BFD_LABELS_DB`, attaches it read-only, searches it before `bfd`, and warns when the pipeline file is newer. Nothing uses it by default yet; `db/BFD.duckdb` and `workflows/build_bfd_db.sh` are unchanged until a fixed pipeline build exists. Checked on the live DB (job 32339225, 24,363 rows): same 43 columns and rows as the old view, except the 8 Antarctolithica genomes, which get Antarctolithicomycotina (SUBPHYLUM NULL 30 to 22).
- config/refs/species_subphylum.tsv, `ref_species_subphylum` and the `species_table` value of `subphylum_source`: a SPECIES-name placeholder for genomes with no CLASS, used after the NCBI value and before the class and phylum tables. First row: `Ascomycota sp. Antarctolithica obscura` (8 genomes) to `Antarctolithicomycotina` (user's request 2026-11-21). On the live DB, genomes without a SUBPHYLUM should fall from 30 to 22 after the next `build_bfd_db.sh` run (not yet run). `build_bfd_db.sh` now requires the new file.
- config/refs: SUBPHYLUM mappings added at the user's request 2026-11-19: eight NCBI classes (Chytridiomycetes, Entomophthoromycetes, Peribolosporomycetes, Novakomycetes, Xylobotryomycetes, Caulochytriomycetes, Olpidiomycetes, Calcarisporiellomycetes), Neocallimastigomycetes (class), Sanchytriomycetes, Monoblepharidomycetes and Hyaloraphidiomycetes (class), and Basidiobolomycota, Blastocladiomycota, Microsporidia (subphylum set to the phylum name) Rozellomycota and Aphelidiomycota (phylum fallback). On the live DB, genomes without a SUBPHYLUM fall from 267 to 30; no genome with an NCBI value changes. Takes effect at the next `build_bfd_db.sh` run.
- `v_genome_labels` carries FungalTraits 1.2 columns (`growth_form`, `growth_form_source`, `primary_lifestyle`,
  `secondary_lifestyle`, `plant_pathogenic_capacity`, `decay_substrate`, `animal_biotrophic_capacity`,
  `in_fungaltraits`, `yeast_origin`) and `SUBPHYLUM` / `subphylum_source`. `sql/bfd_traits.sql` loads the reference
  files in `config/refs/` into three `ref_*` tables and defines `v_genome_traits`; `workflows/build_bfd_db.sh` applies it
  and checks the result before publishing. `bin/check_traits_vs_atlas.py` compares the DB with the atlas parquet.
  `SUBPHYLUM` is NULL when the class has no entry in `config/refs/class_subphylum.tsv` (rows marked `proposed` are unchecked).
- `pyproject.toml` (PEP 621) describing the project and dependencies: core
  (`duckdb`, `numpy`, `pandas`, `scipy`, `PyYAML`) plus `phylogeny`
  (`biopython`) and `alt-ordinations` (`scikit-learn`, `umap-learn`) extras.
- Initial release of the generalized `BFD_dashboard` framework (port of the
  one-off `Rhodotorula/Novelty/dashboard`).

### Feature sets (7)
- `assembly` — assembly/annotation stats (size, N50, GC%, gene count, BUSCO), z-scored, Euclidean + PCoA.
- `codon` — codon usage frequencies (64), chi-square + correspondence analysis.
- `amino_acid` — amino-acid usage (20), chi-square + correspondence analysis.
- `pfam` / `cazy` / `merops` — domain families, gene-count-normalized per 1000 genes, Bray-Curtis + PCoA.
- `localization` — signalp/tmhmm/targetp/predgpi/wolfpsort fractions, z-scored, Euclidean + PCoA.

### Subsetting
- Taxonomy filters: `--phylum`, `--class`, `--order`, `--family`, `--genus`, `--species` (substring match).
- `--only-with-functional` — keep only genomes present in a functional-annotation table.
- `--max-per-species N` — best-assembled strain per species (gene count, then N50).
- `--max-genomes N` with `--seed` — reproducible random cap.
- Empty subsets abort with an explanation; single-genome subsets warn.

### Configuration
- Single `DashboardConfig`/`SubsetConfig` surface; YAML configs layered under CLI flags.
- `config/BFD_dashboard.yaml`, `config/rhodotorula_full.yaml`, `config/rhodotorula_reps.yaml`, `config/full_db_genus_subset.yaml`.

### Analysis
- PCoA / CA / NMDS / t-SNE / UMAP ordinations, PERMANOVA by genus and per-species
  (auto-skipped above 40 genomes).
- Per-species differential domain-family tables (present-only / absent-only / fold-change) with example-protein drill-down.
- Optional phylogeny panel: phyling FastTree tip order, pruned to subset, outgroup-rooted.

### Robustness against large/sparse databases
- Empty or all-zero feature matrices are skipped with a per-feature-set notice instead of crashing.
- Genome stats fall back from `v_species_summary` to `asm_stats` so a sparse summary doesn't blank the overview.
- MEROPS families file optional; skipped cleanly when unconfigured.
- `lib/ordination.py` `_pseudo_f` returns `nan` for single-group PERMANOVA instead of raising.

## [0.1.0] - 2026-08-10

### Added
- Initial commit of the framework.
