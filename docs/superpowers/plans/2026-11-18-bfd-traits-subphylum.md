# FungalTraits and subphylum in the dashboard DB: Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Put FungalTraits ecology columns and a subphylum column into the dumped `db/BFD.duckdb` and into `v_genome_labels`.

**Architecture:** Reference files in `config/refs/` are loaded into three `ref_*` tables at build time by a new `sql/bfd_traits.sql`. That file also defines `v_genome_traits` (one row per genome). `v_genome_labels` left-joins `v_genome_traits`. `workflows/build_bfd_db.sh` applies `bfd_traits.sql` before `bfd_extra_views.sql` and verifies the result before publishing.

**Tech Stack:** DuckDB SQL (CLI via `module load duckdb`, Python `duckdb` under `/usr/bin/python3.12`), bash, pytest.

**Spec:** `docs/superpowers/specs/2026-11-18-bfd-traits-subphylum-design.md`. One change from the spec: the trait join lives in a separate `v_genome_traits` view (spec "Migration" step 3 offers this) so `v_genome_labels` is edited, not copied.

## Global Constraints

- Work only in worktree `.claude/worktrees/traits-import` (branch `traits-import`, from `37826f1`). Do not edit `Fungi_BFD`, `bin/summary/*`, `docs/quarto_summary/*`, `_quarto.yml`.
- Do not write `db/BFD.duckdb` in the main checkout. Test builds use a different `OUT`.
- Do not rewrite `Fungi_BFD_runs/analysis/summary_tables/genome_all.parquet` (SLURM job 32173269 reads it).
- Column names and rules equal `bin/atlas/build_labels.py`: `growth_form`, `growth_form_source`, `primary_lifestyle`, `secondary_lifestyle`, `plant_pathogenic_capacity`, `decay_substrate`, `animal_biotrophic_capacity`, `in_fungaltraits`, `yeast_origin`. New: `SUBPHYLUM`, `subphylum_source`.
- No phylum fallback in the DB: unmapped CLASS gives NULL `SUBPHYLUM`.
- `/usr/bin/python3.12` for Python. No `BASH_SOURCE[0]` in scripts (SLURM). Use absolute paths or environment variables.
- Commit messages end with `Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>`.
- Cite FungalTraits as Polme et al. 2020, Fungal Diversity 105:1-16, marked "verify".

## Review Focus

- FungalTraits has 10,771 rows for 10,765 genera: duplicate genera must not multiply `v_genome_labels` rows, and the kept row must be deterministic (first by `jrk_template`). Pinned in Task 2.
- `species.GENUS` NULL: the genome keeps one row with NULL traits. Pinned in Task 2.
- Override file row with empty `growth_form`: ignored (`nullif`), FungalTraits value kept. Pinned in Task 2.
- `species.SUBPHYLUM` already filled (after a future source fix): it wins over the class table, `subphylum_source = 'species'`. Pinned in Task 2.
- CLASS NULL or absent from the table: `SUBPHYLUM` NULL, `subphylum_source` NULL. Pinned in Task 2.
- Build run while the old DB exists: failure leaves the old DB in place. Already true of the script; Task 4 keeps the checks before the `mv`.

---

### Task 1: Reference files

**Files:**
- Create: `config/refs/fungaltraits_1.2.utf8.csv`, `config/refs/growth_form_overrides.tsv`, `config/refs/class_subphylum.tsv`, `config/refs/README.md`
- Test: `tests/test_traits_refs.py`

**Interfaces:**
- Produces: `class_subphylum.tsv` with header `CLASS\tSUBPHYLUM\tsource\tnote`. `source` is `taxonomy_py` (copied from `bin/summary/taxonomy.py`), `atlas` (from the lists in `bin/atlas/build_labels.py`), or `proposed` (new; user must check). One row per CLASS.

- [ ] **Step 1: Write the failing test** in `tests/test_traits_refs.py`: `test_ft_csv_is_utf8_with_genus_column` (reads with `encoding="utf-8"`, header contains `GENUS` and `primary_lifestyle`, at least 10,000 rows); `test_class_subphylum_has_header_and_unique_classes`; `test_class_subphylum_sources_are_known` (values in `{taxonomy_py, atlas, proposed}`); `test_data_spellings_are_mapped` (`Pneumocystomycetes`, `Candelariomycetes`, `Spiculogloeomycetes`, `Moniliellomycetes` are present: the data spells the first one without "id", unlike `taxonomy.py`).
- [ ] **Step 2: Run it and confirm it fails.** `cd <worktree> && /usr/bin/python3.12 -m pytest tests/test_traits_refs.py -q` fails: files missing.
- [ ] **Step 3: Create the files.**
  - FungalTraits: `iconv -f latin1 -t utf-8 /bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/refs/lifestyle/FungalTraits_1.2.csv > config/refs/fungaltraits_1.2.utf8.csv`; then `cmp` it with `/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/atlas_data/labels/fungaltraits_1.2.utf8.csv` (must be identical).
  - Overrides: copy `/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/atlas_data/labels/growth_form_overrides.tsv` unchanged.
  - `class_subphylum.tsv`: every row of `CLASS_TO_SUBPHYLUM` in `bin/summary/taxonomy.py` with `source=taxonomy_py`, except spelling `Pneumocystidomycetes` is written `Pneumocystomycetes` (the data's spelling; note it); add the classes named in the `SACCH/TAPHRINA/PUCCINIO_YEAST/USTILAGO/PEZIZO` lists of `build_labels.py` that are missing, with `source=atlas`.
  - Then list classes in the real DB that are still unmapped: `duckdb -readonly ../../../db/BFD.duckdb -c "SELECT PHYLUM, CLASS, count(*) FROM species GROUP BY ALL ORDER BY 3 DESC"` (path is the main checkout's DB) and add `proposed` rows only where the subphylum is known (for example Mucoromycetes and Umbelopsidomycetes to Mucoromycotina, Zoopagomycetes to Zoopagomycotina, Wallemiomycetes to Wallemiomycotina). Leave out classes you cannot place; say so in `note`. Do not guess.
  - `README.md`: source, date, citation, and that `proposed` rows are unchecked.
- [ ] **Step 4: Run the test and confirm it passes.**
- [ ] **Step 5: Commit** `config/refs tests/test_traits_refs.py` with message `refs: FungalTraits 1.2, growth form overrides, class-to-subphylum table`.

---

### Task 2: `sql/bfd_traits.sql`

**Files:**
- Create: `sql/bfd_traits.sql`, `tests/test_traits_sql.py`

**Interfaces:**
- Consumes: Task 1 files; DuckDB view/table `species` (columns `LOCUSTAG, GENUS, CLASS, SUBPHYLUM`).
- Produces: tables `ref_fungaltraits`, `ref_growth_form_override`, `ref_class_subphylum`; view `v_genome_traits(LOCUSTAG, growth_form, growth_form_source, primary_lifestyle, secondary_lifestyle, plant_pathogenic_capacity, decay_substrate, animal_biotrophic_capacity, in_fungaltraits, yeast_origin, SUBPHYLUM, subphylum_source)`. Line 1 of the file is `-- REF_DIR=/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/config/refs`; every file path in the SQL is `<REF_DIR>/<name>` using that literal so the build can substitute it.

- [ ] **Step 1: Write the failing tests** in `tests/test_traits_sql.py`. Helper `load(tmp_path, species_rows, ft_rows, ovr_rows, cls_rows) -> duckdb connection`: writes tiny CSV/TSV files into `tmp_path`, creates an in-memory `species` table, replaces the `REF_DIR` literal in `sql/bfd_traits.sql` with `tmp_path`, runs it. Tests:
  - `test_duplicate_genus_keeps_first_by_jrk_template`: two FT rows for `Mucor` with `jrk_template` `jrk2`, `jrk1` and different `primary_lifestyle`; view has one row per genome and the `jrk1` value.
  - `test_override_wins_and_sets_source`: `growth_form_source == 'override'`; an override with empty `growth_form` leaves `'fungaltraits'`.
  - `test_unmatched_genus_and_null_genus`: `in_fungaltraits` false, traits NULL, one row each.
  - `test_subphylum_from_class_table`: `SUBPHYLUM == 'Pezizomycotina'`, `subphylum_source == 'class_table'`; unmapped class and NULL class give NULL and NULL.
  - `test_species_subphylum_wins`: filled `species.SUBPHYLUM` is kept, `subphylum_source == 'species'`.
  - `test_yeast_origin_matches_atlas_rules`: growth_form `yeast` with CLASS `Saccharomycetes` gives `Saccharomycotina`; `dimorphic_yeast` with `Eurotiomycetes` gives `Pezizomycotina-dimorphic`; `filamentous_mycelium` gives `filamentous`; other gives `other`.
- [ ] **Step 2: Run and confirm failure** (`sql/bfd_traits.sql` missing).
- [ ] **Step 3: Write `sql/bfd_traits.sql`.** Create the three tables with `CREATE OR REPLACE TABLE ... AS SELECT ... FROM read_csv(..., header=true, all_varchar=true)` (overrides: `delim='\t'`, `nullif(growth_form,'')`). `ref_fungaltraits` has one row per `GENUS`: `row_number() OVER (PARTITION BY GENUS ORDER BY jrk_template) = 1`, with the columns and renames of the `ft` table in `build_labels.py` lines 85-92. `v_genome_traits` follows `build_labels.py` lines 96-131 for `growth_form`, `growth_form_source`, `in_fungaltraits` and `yeast_origin` (copy the CASE and the class lists from lines 21-30 into SQL), and `SUBPHYLUM = coalesce(s.SUBPHYLUM, c.SUBPHYLUM)` with `subphylum_source` `'species'`, `'class_table'` or NULL. Add `COMMENT ON` for each table and the view.
- [ ] **Step 4: Run the tests and confirm they pass.**
- [ ] **Step 5: Commit** `sql/bfd_traits.sql tests/test_traits_sql.py` with message `sql: ref tables and v_genome_traits for FungalTraits and subphylum`.

---

### Task 3: `v_genome_labels` joins the traits

**Files:**
- Modify: `sql/bfd_extra_views.sql` (`v_genome_labels`, about lines 39-84: select list, `COMMENT`, the "Not included" note)
- Test: `tests/test_traits_sql.py`

**Interfaces:**
- Consumes: `v_genome_traits` from Task 2.
- Produces: `v_genome_labels` with the eleven new columns from the Global Constraints (nine trait columns, `SUBPHYLUM`, `subphylum_source`), one row per `species` row, `LOCUSTAG` unique.

- [ ] **Step 1: Write the failing test** `test_v_genome_labels_one_row_per_genome_with_traits`: build synthetic `species`, `samples`, `asm_stats`, `busco_genome` tables (columns used by the view: `species(LOCUSTAG, ASMID, PHYLUM, CLASS, "ORDER", FAMILY, GENUS, SPECIES, STRAIN, NCBI_TAXONID, SUBPHYLUM)`, `samples(LOCUSTAG, TRANSL_TABLE)`, `asm_stats(ASMID, annotation_status, TOTAL_LENGTH, contig_count, N50_bp, GC_PERCENT, masked_pct)`, `busco_genome(ASMID, complete_pct, single_pct, duplicated_pct, fragmented_pct, missing_pct)`), apply `bfd_traits.sql` as in Task 2, then execute only the `v_genome_labels` statement (slice `sql/bfd_extra_views.sql` from `CREATE OR REPLACE VIEW v_genome_labels AS` to the first `;`; the other views need tables the fixture does not have). Assert `count(*)` equals the `species` count, `count(DISTINCT LOCUSTAG)` equals it, and the trait columns are present with the expected values for one genome. Add a `Mucor` duplicate in the FT fixture so a join that multiplies rows fails the count.
- [ ] **Step 2: Run and confirm it fails.**
- [ ] **Step 3: Edit `v_genome_labels`:** add `LEFT JOIN v_genome_traits t ON s.LOCUSTAG = t.LOCUSTAG` in `base` and select `t.*` minus `LOCUSTAG` (`t.* EXCLUDE (LOCUSTAG)`). Update the comment block and the `COMMENT ON VIEW` text: the FungalTraits and subphylum columns are included, `SUBPHYLUM` is NULL when unmapped, `in_fungaltraits` false means the genus is not in FungalTraits 1.2.
- [ ] **Step 4: Run all tests in `tests/` and confirm they pass.**
- [ ] **Step 5: Commit** with message `sql: v_genome_labels carries FungalTraits columns and SUBPHYLUM`.

---

### Task 4: Build script wiring and checks

**Files:**
- Modify: `workflows/build_bfd_db.sh`
- Create: `bin/check_traits_vs_atlas.py`

**Interfaces:**
- Consumes: `sql/bfd_traits.sql` with `-- REF_DIR=` first line.
- Produces: `bin/check_traits_vs_atlas.py --db PATH --atlas PARQUET` prints per-column mismatch counts for `growth_form`, `primary_lifestyle`, `SUBPHYLUM` null count, and unmapped classes by genome count. Exit 1 if `growth_form` or `primary_lifestyle` differ for any genome whose genus has no duplicate FungalTraits row; exit 0 otherwise. Reads compressed or plain parquet through DuckDB.

- [ ] **Step 1: Write the failing test** `test_check_script_flags_a_mismatch` in `tests/test_traits_sql.py`: tiny DB file and tiny atlas parquet with one differing `primary_lifestyle`; run the script through `subprocess` and assert exit code 1 and the genome id in stdout. Add `test_check_script_passes_when_equal`.
- [ ] **Step 2: Run and confirm it fails.**
- [ ] **Step 3: Implement `bin/check_traits_vs_atlas.py`.**
- [ ] **Step 4: Edit `build_bfd_db.sh`.** Add `REF_DIR="${REF_DIR:-$REPO/config/refs}"` next to the other variables. Before the `bfd_extra_views.sql` loop iteration, apply `bfd_traits.sql`: substitute the line-1 `REF_DIR` literal with `$REF_DIR`, as the loop does for `TABLES`, and run it with `duckdb "$STAGE/BFD.duckdb"`. Order: `bfd_compat_views.sql`, `bfd_traits.sql`, `bfd_extra_views.sql`. In the Python verify block: add `v_genome_traits` to the helper-view list; check the three `ref_*` tables are non-empty; check `count(*) FROM v_genome_labels` equals `count(*) FROM species` and `count(DISTINCT LOCUSTAG)` equals it; print the number of genomes with NULL `SUBPHYLUM`. All checks append to `problems` so a failure exits before the `mv`. Run `bash -n workflows/build_bfd_db.sh`.
- [ ] **Step 5: Run the tests and `bash -n`; confirm they pass.**
- [ ] **Step 6: Commit** with message `build: apply bfd_traits.sql and verify trait columns before publish`.

---

### Task 5: Real-data check, documentation, handoff

**Files:**
- Modify: `CHANGELOG.md`, `docs/BFD_VIEW_REVIEW.md` (the `v_genome_labels` and "Not included" statements), `docs/superpowers/specs/2026-11-18-bfd-traits-subphylum-design.md` (record the `v_genome_traits` change)

- [ ] **Step 1: Check on a copy of the real DB.** `cp ../../../db/BFD.duckdb $SCRATCH/BFD_traits_copy.duckdb`; apply `bfd_traits.sql` with `REF_DIR` substituted to this worktree's `config/refs` (use `sed` as in the build script) and the new `v_genome_labels` statement; then run `/usr/bin/python3.12 bin/check_traits_vs_atlas.py --db $SCRATCH/BFD_traits_copy.duckdb --atlas /bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/atlas_data/labels/genome_labels.parquet`. Expected: exit 0; `count(*)` of `v_genome_labels` is 23,680. Report the exact mismatch and unmapped-class numbers. Do not call it equal unless the script says so.
- [ ] **Step 2: Full build to a separate path.** Confirm no `MERGE_*` or Nextflow run is writing `Fungi_BFD_runs/tables/` (`squeue -u $USER`, `pgrep -af nextflow`), then `OUT=<worktree>/db_test/BFD.duckdb sbatch workflows/build_bfd_db.sh`. Report the verify output. Do not set `OUT` to `db/BFD.duckdb`.
- [ ] **Step 3: Docs.** Add a `CHANGELOG.md` entry; update the two docs; update the spec line about `v_genome_labels` replacement.
- [ ] **Step 4: Run all tests; commit** with message `docs: changelog and view review for traits and subphylum`.
- [ ] **Step 5: Push the branch** over `https://github.com/stajichlab/BFD_dashboard.git` and open a PR whose body lists the check numbers from Steps 1 and 2 and the `proposed` subphylum rows for the user to check. Do not merge.
