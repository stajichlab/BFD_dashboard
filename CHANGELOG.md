# Changelog

All notable changes to this project are documented in this file.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

The version number is bumped when a feature or fix is added. Until the project
reaches a stable 1.0, minor version bumps mean "something was added or changed".

## [Unreleased]

### Added
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
