# Ecology in the summary ordinations: Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A new step, `bin/summary/summary_ecology.py`, that draws the site's genus-level ordinations coloured by FungalTraits lifestyle and growth form, tests lifestyle on the whole tree and inside every qualifying clade, and a page `17_ecology.qmd` that shows the results.

**Architecture:** One new module with small pure functions (permutation test, clade selection, plotting) and a `main()` that wires them. It imports the PCoA and matrix code from `summary_ordination.py` and `lib/ordination.py` and the `Out` class from `make_figures.py`, and changes none of them. Traits come from `v_genome_labels` in the DB. Outputs go to their own directory.

**Tech Stack:** Python 3.12 (`/usr/bin/python3.12`), pandas, numpy, matplotlib, duckdb, pytest. Bash for the SLURM script. Quarto page.

**Spec:** `docs/superpowers/specs/2026-11-18-bfd-ecology-ordination-design.md`, including its section "Calls made after the builder check". The spec governs wherever this plan is silent.

## Global Constraints

- Work in worktree `.claude/worktrees/ecology-ordination` (branch `ecology-ordination`, from `133eedb`).
- All new files only. Do not edit: `bin/summary/make_figures.py`, `summary_ordination.py`, `figure_index.py`, `lib/ordination.py`, `docs/quarto_summary/_quarto.yml`, `06_data_notes.qmd`, `figures_release.json`, `tf_*`, `genome_all.parquet`. Read `analysis/summary_tables/*.parquet` only.
- No release upload. No push or PR without the user's go-ahead.
- Defaults: `--min-genera 30`, `--min-group 5`, `--levels phylum class order`, 999 permutations, seed 1 (`so.N_PERM`, `so.SEED`).
- Clade rule, applied to the ordination reps from one `pick_complete_reps(genome_all, mats)` call with all matrices present, labeled genera only: at least `min_genera` labeled genera and at least 2 lifestyles with at least `min_group` genera each.
- Stratum for the stratified test: ORDER for phylum and class clades, FAMILY for order clades. NULL stratum: left out of the stratified test only and counted. Stratified result is NA with a reason when fewer than 2 strata have at least 2 genera each.
- Test genera whose lifestyle has fewer than `min_group` genera in the tested set are dropped from the test and counted.
- Benjamini-Hochberg q-values are computed per (level, matrix, test).
- Output names: figures `E01_<matrix>_pcoa_lifestyle`, `E02_<matrix>_pcoa_growthform`, `E03_<level>_<slug>_<matrix>_pcoa_lifestyle`; tables `ecology_traits_coverage`, `ecology_permanova_global`, `ecology_permanova_clades`, `ecology_clades_skipped` (all `.tsv`); topic `"ecology"`; figures under `OUT/figures/ecology`, tables under `OUT/tables/ecology`.
- Python: `/usr/bin/python3.12 -m pytest bin/summary/test_summary_ecology.py -q`. No `BASH_SOURCE[0]` in scripts (SLURM). Commit messages end with the two lines `Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>` and `Claude-Session: https://claude.ai/code/session_01CoSBcNwNQyDtYjjwF`.

## Review Focus

- A clade whose labeled genera fall below the threshold only after NULL labels are dropped goes to the skipped table with a reason; it never raises. Pinned in Task 2.
- A clade (or the whole set) where all labeled genera share one lifestyle: skipped or NA with a reason, not an error. Pinned in Tasks 1 and 4.
- Genera with NULL ORDER or FAMILY: counted and left out of the stratified test, kept in the unrestricted one. Pinned in Task 1.
- 28 lifestyle values: pooling keeps the legend inside the canvas and the saved PNG unclipped. Pinned in Task 3.
- Clade names with spaces, slashes or non-ASCII characters, and the same name at two levels: file names stay safe and unique. Pinned in Task 4.
- A matrix parquet file missing: that matrix is skipped and the others run. An old DB without trait columns: exit 1 with a message that names the missing column. Pinned in Task 5.

---

### Task 1: Statistics helpers

**Files:**
- Create: `bin/summary/summary_ecology.py`, `bin/summary/test_summary_ecology.py`

**Interfaces:**
- Consumes: `summary_ordination` as `so` (`so._permanova`, `so._r2_adj`, `so.lib_ordination._pseudo_f(dist: np.ndarray, labels: np.ndarray) -> float`, `so.N_PERM`, `so.SEED`).
- Produces in `summary_ecology.py`:
  - `permute_within(labels: np.ndarray, strata: np.ndarray, rng: np.random.Generator) -> np.ndarray` (labels permuted inside each stratum only).
  - `stratified_permanova(dist: pd.DataFrame, groups: pd.Series, strata: pd.Series | None = None, n_permutations: int = so.N_PERM, seed: int = so.SEED) -> dict` with keys `F, R2, R2_adj, p, n, n_groups, n_strata_used, n_dropped_no_stratum, n_permutations_run, caveat`. `strata=None` permutes all labels. p = (1 + count of permuted F >= observed F) / (1 + n_permutations).
  - `bh_adjust(p: pd.Series) -> pd.Series` (NaN stays NaN).
  - `pool_groups(labels: pd.Series, min_group: int, other: str = "other") -> pd.Series` (groups under `min_group` become `f"{other} ({k} groups)"`; NULL stays NULL).

- [ ] **Step 1: Write the failing tests** in `test_summary_ecology.py` (import `summary_ecology as se`, `summary_ordination as so`; build data with a local `two_cluster_dist()` helper of 40 points):
  - `test_unstratified_matches_existing_permanova`: `se.stratified_permanova(dist, groups)` gives `F` and `R2` equal to `so._permanova(dist, groups)` (`pytest.approx`, rel 1e-9), and `p < 0.01` on separated clusters.
  - `test_permute_within_keeps_each_stratum_multiset`: over 50 draws, sorted labels inside every stratum are unchanged.
  - `test_confounded_labels_give_p_one_when_stratified`: labels constant within every stratum; unrestricted `p < 0.01`, stratified `p == 1.0`.
  - `test_fewer_than_two_informative_strata_is_na`: `p` is NaN, `caveat` mentions "strata".
  - `test_null_stratum_rows_counted_and_kept_out`: `n_dropped_no_stratum == 3` for 3 NULL strata; the unrestricted call keeps them in `n`.
  - `test_single_group_is_na_not_error`: one lifestyle only gives NaN and a caveat.
  - `test_bh_adjust_known_values`: `[0.01, 0.02, 0.03]` gives `[0.03, 0.03, 0.03]`; NaN preserved.
  - `test_pool_groups`: groups of size 4 and 2 with `min_group=5` become `"other (2 groups)"`; NULL stays NULL.
- [ ] **Step 2: Run and confirm failure** (`ModuleNotFoundError` for `summary_ecology`).
- [ ] **Step 3: Implement the four functions.** R2 and R2_adj use the same formulas as `so._permanova`; the observed and permuted F come from `lib_ordination._pseudo_f`. Strata with fewer than 2 rows do not move.
- [ ] **Step 4: Run the tests; confirm all pass.**
- [ ] **Step 5: Commit** `summary_ecology.py` and the test file: `ecology: stratified PERMANOVA, BH, group pooling`.

---

### Task 2: Traits and clade selection

**Files:**
- Modify: `bin/summary/summary_ecology.py`, `bin/summary/test_summary_ecology.py`

**Interfaces:**
- Produces:
  - `load_traits(db_path: str) -> pd.DataFrame` with columns `species_prefix, primary_lifestyle, growth_form`, read-only from `v_genome_labels` (`LOCUSTAG AS species_prefix`). If a column is missing it raises `RuntimeError` whose text names the missing column and says to run `workflows/build_bfd_db.sh`.
  - `@dataclass Clade: level: str; name: str; ids: list[str]`
  - `select_clades(reps: pd.DataFrame, labels: pd.Series, levels: tuple[str, ...] = ("phylum", "class", "order"), min_genera: int = 30, min_group: int = 5) -> tuple[list[Clade], pd.DataFrame]`. `reps` has `species_prefix, PHYLUM, CLASS, ORDER`; `labels` is indexed by `species_prefix`. The second value is the skipped table with columns `level, clade, n_labeled, n_groups_ge_min, reason`. Clades with a NULL name are not listed.
  - `LEVEL_COLUMN = {"phylum": "PHYLUM", "class": "CLASS", "order": "ORDER"}` and `STRATUM_COLUMN = {"phylum": "ORDER", "class": "ORDER", "order": "FAMILY"}`.

- [ ] **Step 1: Write the failing tests:** `test_select_clades_threshold_and_reasons` (hand-built reps: clade A passes; B has 29 labeled genera, reason contains "genera"; C has 40 genera but one lifestyle, reason contains "lifestyles"; D has 31 genera of which 2 labeled, so the reason counts labeled genera only); `test_select_clades_drops_null_labels_before_counting`; `test_load_traits_from_duckdb` (tiny DB file with a `v_genome_labels` table); `test_load_traits_names_missing_column` (table without `primary_lifestyle`: the error names it).
- [ ] **Step 2: Run and confirm failure.**
- [ ] **Step 3: Implement** the two functions, the dataclass and the two maps.
- [ ] **Step 4: Run the tests; confirm all pass.**
- [ ] **Step 5: Commit:** `ecology: trait loading and clade selection`.

---

### Task 3: Global figures and tables

**Files:**
- Modify: `bin/summary/summary_ecology.py`, `bin/summary/test_summary_ecology.py`

**Interfaces:**
- Consumes: Task 1 and 2 functions; `so._feature_matrix(name, table, reps) -> (matrix, n_missing, n_nan)`, `so.pick_complete_reps`, `so.METRIC`, `so.N_AXES`, `so.lib_ordination.distance_matrix` and `pcoa`; `mf.Out(root, "ecology")` with `.table(df, name)`, `.save(fig, name, topic)`, `.made`; `figs_common.safe`.
- Produces:
  - `lifestyle_palette(groups: list[str]) -> dict[str, str]`: one colour per group, all distinct; the reserved greys `NO_LABEL = "#d0d0d0"` and `OTHER_GREY = "#9a9a9a"` are never used for a real group.
  - `scatter_labels(coords: pd.DataFrame, labels: pd.Series, palette: dict[str, str], pct: list[float], title: str) -> matplotlib.figure.Figure`: legend right of the axes, 2 columns, 6 pt, each entry `"<group> (<n>)"`, unlabeled drawn in `NO_LABEL` and listed as "no label (n)".
  - `run_global(reps: pd.DataFrame, mats: dict[str, pd.DataFrame], traits: pd.DataFrame, out: Path, min_group: int = 5) -> dict`: writes `E01_*` and `E02_*` figures, `ecology_traits_coverage.tsv` (columns `trait, n_genera, n_labeled, n_unlabeled`), `ecology_permanova_global.tsv` (columns `matrix, trait, test, n_genera, n_groups, n_dropped_rare_group, n_dropped_no_stratum, n_strata_used, pseudo_F, R2, R2_adj, p_value, caveat`; `test` is `unrestricted` or `within_order`). Returns per-matrix coordinates `{matrix: DataFrame}`.

- [ ] **Step 1: Write the failing tests** with a synthetic fixture (two phyla, 40 genera, lifestyles and orders assigned by hand): `test_global_coords_match_existing_ordination_up_to_sign` (run `so.run_ordination` and `se.run_global` on the same fixture into two temp dirs; for PCo1 and PCo2, `np.allclose(abs(a), abs(b), atol=1e-9)`); `test_coverage_reconciles` (`n_labeled + n_unlabeled == n_genera` in the table and in the E01 legend counts); `test_global_permanova_has_both_tests`; `test_scatter_legend_outside_and_png_not_clipped` (as the PR #19 tests: legend `x0 >= ax.x1`; saved PNG has white margin on the right) with 28 lifestyle values pooled by `min_group`; `test_palette_distinct_and_reserved_greys_unused`; `test_missing_matrix_is_skipped` (only `codon` in `mats`).
- [ ] **Step 2: Run and confirm failure.**
- [ ] **Step 3: Implement** the three functions. The distance, PCoA and feature-matrix calls are those of `run_ordination`.
- [ ] **Step 4: Run the tests; confirm all pass.**
- [ ] **Step 5: Commit:** `ecology: global ordinations coloured by lifestyle and growth form`.

---

### Task 4: Clade zoom-ins

**Files:**
- Modify: `bin/summary/summary_ecology.py`, `bin/summary/test_summary_ecology.py`

**Interfaces:**
- Consumes: `Clade`, `select_clades`, `stratified_permanova`, `bh_adjust`, `pool_groups`, `scatter_labels`, `lifestyle_palette`.
- Produces:
  - `clade_slug(level: str, name: str) -> str`: lower-case, non-alphanumeric runs become `_`, prefixed with the level (`class_agaricomycetes`); distinct for the same name at two levels.
  - `run_clades(reps: pd.DataFrame, mats: dict[str, pd.DataFrame], traits: pd.DataFrame, clades: list[Clade], skipped: pd.DataFrame, out: Path, min_group: int = 5) -> pd.DataFrame`: for each clade and matrix, recompute the feature matrix (structure z-scoring on the clade's genera only), distance and PCoA on the clade's reps (no re-picking of reps); save `E03_<slug>_<matrix>_pcoa_lifestyle`; run the unrestricted and the stratified tests (`test` is `unrestricted` or `stratified`). Writes `ecology_permanova_clades.tsv` (columns `level, clade, matrix, test, n_genera, n_labeled, n_groups, n_dropped_rare_group, n_dropped_no_stratum, n_strata_used, pseudo_F, R2, R2_adj, p_value, q_BH, caveat`) and `ecology_clades_skipped.tsv`; returns the clade table.

- [ ] **Step 1: Write the failing tests:** `test_clade_slug_safe_and_unique` (`"Nectar/tap sp."`, `"Ünï"`, and the same name at two levels); `test_clade_with_separated_lifestyles_is_significant` (unrestricted `p < 0.01`, `q_BH >= p`); `test_confounded_clade_stratified_p_is_one` (lifestyle follows order inside a class: unrestricted `p < 0.05`, stratified `p == 1.0`); `test_q_bh_is_per_level_matrix_test` (two levels, q computed within each); `test_below_threshold_clade_in_skipped_table`; `test_e03_files_exist_with_expected_names`; `test_structure_scaling_uses_clade_only` (the z-scores of the clade's structure matrix have mean 0 over the clade, not over all genera).
- [ ] **Step 2: Run and confirm failure.**
- [ ] **Step 3: Implement** `clade_slug` and `run_clades`. A clade or matrix with fewer than 3 genera is skipped with a reason in the skipped table.
- [ ] **Step 4: Run the tests; confirm all pass.**
- [ ] **Step 5: Commit:** `ecology: per-clade ordinations and lifestyle tests`.

---

### Task 5: Command line, SLURM script, page

**Files:**
- Modify: `bin/summary/summary_ecology.py`, `bin/summary/test_summary_ecology.py`
- Create: `workflows/summary_ecology.sh`, `docs/quarto_summary/17_ecology.qmd`

**Interfaces:**
- Consumes: Tasks 1 to 4.
- Produces: `main(argv: list[str] | None = None) -> int`. Arguments `--tables DIR --db FILE --out DIR [--min-genera 30] [--min-group 5] [--levels phylum class order]`. It reads `genome_all`, `codon_genome`, `aa_genome` and `gene_structure_genome` parquet files from `--tables` (a missing matrix file is skipped with a printed note; no matrices at all gives exit 1), calls `so.pick_complete_reps(genome_all, mats)` once, merges traits on `species_prefix`, then `run_global`, `select_clades`, `run_clades`. It prints the number of qualifying clades and outputs. Exit 0 on success; exit 1 on any failure, with the message on stderr. Constants `GLOBAL_TABLES`, `CLADE_TABLES` list the table names.
- `workflows/summary_ecology.sh`: SLURM (`-p epyc -c 4 --mem 32G --time 3:00:00 -J summary_ecology`), variables `ROOT` (default `${SLURM_SUBMIT_DIR:?}`), `DATA` (default `/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/summary_tables`), `DB` (default `/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/db/BFD.duckdb`), `OUT` (default `/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/summary_ecology`), `PY` (as `workflows/summary_figures.sh`), `: "${SCRATCH:?}"`, `MPLCONFIGDIR="$SCRATCH/mplconfig"`; runs `$PY bin/summary/summary_ecology.py --tables "$DATA" --db "$DB" --out "$OUT"`.
- `17_ecology.qmd`: same helper import block as `10_ordination.qmd`; text states the label is genus-level from FungalTraits 1.2, the coverage numbers come from `ecology_traits_coverage`, taxonomy confounds lifestyle, and what the two tests mean; sections for the global figures and `ecology_permanova_global`, and the clade table with `ecology_clades_skipped`. How the clade figures are laid out is left to the builder.

- [ ] **Step 1: Write the failing tests:** `test_main_end_to_end_on_synthetic_inputs` (parquet files and a tiny DB file in `tmp_path`; exit 0; every name in `GLOBAL_TABLES` and `CLADE_TABLES` exists under `OUT/tables/ecology`; at least one `E01_` and one `E03_` figure); `test_main_missing_matrix_skips_it`; `test_main_old_db_exits_1_and_names_column` (DB without `primary_lifestyle`: exit 1, stderr names the column); `test_workflow_script_has_no_bash_source_and_parses` (`bash -n`, text lacks `BASH_SOURCE`, contains `SCRATCH` and the four default paths above); `test_page_mentions_every_output_table` (the qmd text contains each name in `GLOBAL_TABLES + CLADE_TABLES`).
- [ ] **Step 2: Run and confirm failure.**
- [ ] **Step 3: Implement** `main`, the script and the page.
- [ ] **Step 4: Run the whole `bin/summary` and `tests` suites in the background** (about 5 minutes): `/usr/bin/python3.12 -m pytest bin/summary tests -q > $SCRATCH/eco_all.log 2>&1`. Expected: all pass. If an existing site test fails because it needs a navbar entry for the new page, rule on the smallest fix that does not edit `_quarto.yml` and write it in the ledger.
- [ ] **Step 5: Run the step on the scratch DB:** `python bin/summary/summary_ecology.py --tables /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/summary_tables --db /bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/db/BFD_traits_dev.duckdb --out $SCRATCH/eco_dev` as a SLURM job (`sbatch` with the new script and `DB`, `OUT` overrides). Expected: exit 0; the printed number of qualifying clades; report the coverage and the clade count. Do not call the result an answer to the ecology question.
- [ ] **Step 6: Commit:** `ecology: command line, SLURM script and page 17`.
