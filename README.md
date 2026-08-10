# BFD Dashboard

A generic, self-contained comparative-genomics dashboard builder for any
`BFD.duckdb` produced by the [BFD nextflow pipeline](https://github.com/stajichlab/BFD) and [Fungi BFD nextflow pipeline](https://github.com/stajichlab/Fungi_BFD). Reads the database and emits a single browsable `explorer.html` (no
server, no CDN) with per-genome overview tiles, ordination + PERMANOVA for a set of
feature matrices, composition bars, and per-species differential domain-family tables.

This is the generalized extraction of the one-off `Rhodotorula/Novelty/dashboard`:
every per-dataset knob (db path, outgroup genus, treefile, title, taxon subset) now
lives in a config file instead of being hardcoded, and any subset of a large database
(e.g. the ~23k-genome full BFD) can be selected by taxonomy and/or subsampled before
the heavy analysis runs.

## Layout

```
BFD_dashboard/
  bin/make_dashboard.py    # CLI: compute payload -> cache/, render -> explorer.html
  lib/
    bfd_data.py            # duckdb -> per-species feature matrices (subset-aware)
    ordination.py          # PCoA / CA / NMDS / t-SNE / UMAP / PERMANOVA
    subset.py              # taxonomy filters + subsampling -> LOCUSTAG list
    figures.py             # CVD-safe palette (self-contained)
    config.py              # DashboardConfig / SubsetConfig dataclasses + YAML load
  config/
    BFD_dashboard.yaml         # example: full 23k-genome db, no subset
    rhodotorula_full.yaml      # focused 19-genome run (validated reference build)
    rhodotorula_reps.yaml      # 1 representative per Rhodotorula species
    full_db_genus_subset.yaml  # genus=Fusarium, capped at 30 genomes (seeded)
  templates/explorer_template.py   # HTML/JS/SVG renderer
  cache/payload.json         # intermediate compute output
  explorer.html              # the dashboard you open in a browser
```

## Quick start

```bash
# Full database, all feature sets, no subset:
python3 bin/make_dashboard.py --config config/BFD_dashboard.yaml

# A subset: genus filter + one representative per species, seeded:
python3 bin/make_dashboard.py --config config/rhodotorula_reps.yaml

# Or express the subset directly on the CLI (overrides config):
python3 bin/make_dashboard.py --db /path/BFD.duckdb \
    --genus Fusarium --max-genomes 30 --seed 7
```

Computation and rendering are two explicit steps so re-theming the HTML doesn't
re-run PERMANOVA:

```bash
python3 bin/make_dashboard.py --compute   # duckdb -> cache/payload.json
python3 bin/make_dashboard.py --render    # cache/payload.json -> explorer.html
```

## Subsetting

Any combination of ranks (case-insensitive substring match on the `species` table) —
`--phylum`, `--class`, `--order`, `--family`, `--genus`, `--species` — plus:

- `--only-with-functional`: keep only genomes present in a functional-annotation
  table (pfam/cazy/merops/signalp/tmhmm/targetp/idp/wolfpsort/predgpi). Essential
  against a large BFD db where those tables are only populated for genomes that have
  finished the pipeline (check `lib/bfd_data.feature_coverage`).
- `--max-per-species N`: at most N strains per species, choosing the best-assembled
  strain (gene count, then N50). `1` = one representative per species.
- `--max-genomes N`: seeded random cap on the total (reproducible for a given `--seed`).

The resolved subset (count + per-filter accounting) is printed to stderr before any
analysis runs. A subset that resolves to 0 genomes aborts with an explanation.

## What is computed

Per feature set (`assembly`, `codon`, `amino_acid`, `pfam`, `cazy`, `merops`,
`localization`):

- distance matrix (Euclidean / Bray-Curtis / chi-square) and a PCoA or correspondence
  analysis of it
- PERMANOVA grouped by GENUS, plus per-species "vs rest" PERMANOVA (auto-skipped above
  ~40 genomes; force with `--per-species-permanova`)
- domain families get a differential table per species: present-only / absent-only /
  top fold-change vs the mean-of-rest, each linking to example `protein_id`s

Pfam/CAZy/MEROPS counts are gene-count-normalized (per 1000 genes) before any
distance step; codon/AA profiles get correspondence analysis; assembly and
localization matrices are z-scored. Zero-handling and normalization choices are
documented once in `lib/bfd_data.py` and enforced consistently.

## The phylogeny panel (optional)

Pass an explicit `--tree <file>` or set `results_dir` in the config to glob for a
phyling FastTree treefile
(`phyling_pep/protein/buildtree/*/fasttree/*.fasttree.support.treefile`). The tree is
pruned to the selected subset and, when `--outgroup-genus` names a GENUS present in
both the tree and the subset, rooted on it. Species ordering across the dashboard
follows the tree's tip order when a tree loads; otherwise it falls back to taxonomic
order (the tree panel then says "no phylogeny available for this build"). A subset
that excludes the outgroup simply falls back to taxonomic order.

## Notes on large databases

The builder:

- skips any feature set whose matrix is empty/all-zero for the chosen subset, with a
  per-feature-set "no data" notice, instead of crashing;
- reports per-table species coverage (`feature_coverage`) so you can see exactly which
  views will be empty;
- generates the caveats text (`notes` in the payload) from the actual subset size and
  genus composition, never a hardcoded "11 genomes" line.

Genome stats in the Overview tiles come from `asm_stats` (available for nearly all
genomes) with `v_species_summary` values preferred where present, so a sparse summary
table doesn't blank out the whole dataset.

## Dependencies

The lab base conda `python3` already has `duckdb`, `numpy`, `scipy`, `pandas`,
`scikit-learn`, `biopython`, `pyyaml`. `umap-learn` is optional (only for
`--alt-ordinations`).
