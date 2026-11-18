# Design: FungalTraits and subphylum in the dashboard DB build

Date: 2026-11-18. Status: approved by the user in conversation (dashboard only now; migration to Fungi_BFD later).

## Intent
FungalTraits ecology (growth form, lifestyle) and subphylum must be available in the dumped
`db/BFD.duckdb` and in `v_genome_labels`. Then SQL users and the Quarto ordinations can colour or group by
ecology without reading files outside the database.

## Facts (measured 2026-11-18)
- `SUBPHYLUM` is NULL in all 23,680 rows of `species` (dashboard DB and `Fungi_BFD_runs/tables/species.parquet`).
- `bin/summary/taxonomy.py` `CLASS_TO_SUBPHYLUM` (written from memory, citation unchecked) covers 22,560 of 23,680 genomes. Unmapped: Mucoromycetes 275, Wallemiomycetes 51, no CLASS 48, Umbelopsidomycetes 14, Pneumocystomycetes 11, Zoopagomycetes 11, others smaller.
- FungalTraits columns exist only in `atlas_data/labels/genome_labels.parquet` (built by `bin/atlas/build_labels.py`). `primary_lifestyle` is filled for 21,642 genomes, `growth_form` for 23,321.
- `v_genome_labels` (`sql/bfd_extra_views.sql`) leaves these columns out. `workflows/build_bfd_db.sh` has no import step.
- FungalTraits 1.2 has 10,771 rows for 10,765 genera, so a genus join must be deduplicated.

## Goals
- Three reference tables inside the DB: `ref_fungaltraits`, `ref_growth_form_override`, `ref_class_subphylum`.
- `v_genome_labels` gains: `growth_form`, `growth_form_source`, `primary_lifestyle`, `secondary_lifestyle`,
  `plant_pathogenic_capacity`, `decay_substrate`, `animal_biotrophic_capacity`, `in_fungaltraits`, `yeast_origin`,
  `SUBPHYLUM`, `subphylum_source`. Names and rules match `build_labels.py`.
- The DB equals the atlas parquet for `growth_form` and `primary_lifestyle` on every genome (acceptance test).
- The build verifies row count and key uniqueness and prints the unmapped-class count.

## Non-goals and anti-goals
- No change to `Fungi_BFD` or its checkout in this work.
- No ordination or figure change here. The ecology ordination work is a separate design on the figure side.
- No phylum fallback inside the DB: unmapped class gives NULL `SUBPHYLUM`. Figure code keeps its own fallback.
- Failure case: a join that multiplies rows, or a build that overwrites `db/BFD.duckdb` before verification passes.

## Approaches considered
1. Materialize reference tables from repo files at build time and join in `v_genome_labels`. **Chosen.** The dump is self-contained.
2. Views that read the CSV by path. Rejected: the dump depends on a path.
3. Read `atlas_data/labels/genome_labels.parquet`. Rejected: couples the DB to atlas outputs that can go stale.

## Design
- Reference files in `config/refs/`: `fungaltraits_1.2.utf8.csv` (UTF-8 copy of the DeltaGain_Fungi file), `growth_form_overrides.tsv` (copy of the atlas file), `class_subphylum.tsv` (seeded from `taxonomy.py`, one row per class, with a `source` column).
- New file `sql/bfd_traits.sql`: creates the `ref_*` tables from the files and defines `v_genome_traits` (one row per genome); `v_genome_labels` in `bfd_extra_views.sql` left-joins it. (Changed from "replaces `v_genome_labels`" during planning so the view is edited, not copied.) Repo file paths are substituted at build time the same way `TABLES` is.
- `workflows/build_bfd_db.sh` applies it after `bfd_extra_views.sql` and extends the verify block: `v_genome_labels` count equals `species` count; no duplicate `LOCUSTAG`; the three ref tables are non-empty; print unmapped-class genome count.
- Test, outside the build: compare `growth_form` and `primary_lifestyle` per genome with the atlas parquet.
- FungalTraits is cited as Polme et al. 2020, Fungal Diversity 105:1-16 (verify before publication).
- The proposed subphylum for each unmapped class is listed in `class_subphylum.tsv` with a `source` of `proposed` and a confidence note. The user checks them.
- Testing uses a copy of the DB or an `OUT` path in scratch. `db/BFD.duckdb` is not touched until the build verifies.

## Migration plan to Fungi_BFD (later, separate PR)
1. Keep `sql/bfd_traits.sql` and `config/refs/*` self-contained: no dashboard paths, one `REF_DIR` substitution.
2. After the join is validated here, offer to `Fungi_BFD` a PR that adds the same files and a `mkview`-style step to `nextflow/bin/build_BFD_duckDB.sh`, so `rebuild_pipeline_duckdb.sh` also gets traits.
3. Decide where the views live: the pipeline builder owns `species`; `v_genome_labels` currently lives only in the dashboard SQL. Move it together with the trait join, or add the traits as a separate `v_genome_traits` view that both DBs carry.
4. Fix `SUBPHYLUM` at the source in `Fungi_BFD` (the `species` table is shared by live runs: work in a worktree, never the shared checkout). Then `ref_class_subphylum` becomes a fallback.
5. After the move, `build_bfd_db.sh` stops applying `bfd_traits.sql`. Check by comparing the two DBs' `v_genome_labels`.
Prerequisite: a decision from the user and from the `bfd-build` session on who owns the `Fungi_BFD` change.

## Left to the builder
- Exact SQL for deduplicating FungalTraits by genus: follow `build_labels.py` lines 85-93.
- Whether `yeast_origin` class lists are copied from `build_labels.py` into SQL or into a reference table.
