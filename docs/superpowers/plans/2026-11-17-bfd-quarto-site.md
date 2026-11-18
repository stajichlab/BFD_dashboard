# BFD Quarto Site Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Publish a static Quarto site at `bfd.fungalgenomes.org` with many genome, gene, codon, amino-acid and gene-structure figures and comparisons by phylum, order and subgroup, plus a first interactive (Observable JS) page.

**Architecture:** Extend the existing summary pipeline (`bin/summary/`). `build_summary_tables.py` writes small per-genome parquet tables from the BFD DuckDB. `make_figures.py` and new figure modules write PNG/SVG and TSV per level. Quarto pages in `docs/quarto_summary/` display them. GitHub Actions renders and deploys to Pages. Phase B adds a JSON export and one OJS page.

**Tech Stack:** Python 3.12 (pandas, numpy, scipy, matplotlib, pyarrow, duckdb, pytest), Quarto 1.7.31 (`module load quarto/1.7.31`), GitHub Actions + Pages.

**Spec:** `docs/superpowers/specs/2026-11-17-bfd-quarto-site-design.md`

## Global Constraints

- Repo `stajichlab/BFD_dashboard`. SSH to GitHub is denied on HPCC. Push over `https://github.com/stajichlab/...`.
- Levels: phylum, subphylum, class, order. Groups smaller than `MIN_N` (10) are pooled as `Other (small groups)`.
- Genomes with no genus or order stay in the data. Label them `Unknown` at that level.
- Colors: one hue per phylum. Subgroups are shades of the parent phylum hue. Colorblind-safe.
- Do not use `targetP` or `infernal_rfam`. Do not use tables built before the 2026-11-16 rebuilds without a check.
- Heavy computation runs in SLURM (`epyc`, 16 CPU, 64 GB). Use `$SCRATCH` for DuckDB temp files (`$SCRATCH:?`). Do not touch files under a live Nextflow run.
- Large data stays out of git. `atlas_data/`, `db/` are ignored. The Pages site must stay under 1 GB.
- Tests run with `python3.12 -m pytest -q bin/summary/` in an env that has pandas, numpy, scipy, matplotlib, pyarrow, pytest. The `atlas_env` pixi env has most of these. Add `pytest` if missing.
- Quarto needs `nbformat nbclient ipykernel pyyaml tabulate` in the Python that `QUARTO_PYTHON` points to.
- Compress large text output (zstd or gzip). Keep any new script able to read `.gz` input.

## Review Focus

- All-genome scatter plots: about 23,000 points must stay readable (hexbin or alpha, rasterized points) and must not silently use only the `is_representative` rows or only genus representatives.
- A genome with NULL `GENUS` or `ORDER`: must appear as `Unknown`, not vanish and not crash group sorting.
- Order level has 195 orders: `MAX_GROUPS` (36) and `MIN_N` pooling must keep every phylum visible and keep legends readable.
- A statistic with zero variance in a group, or a group with one genome: ANOVA and regression must return NaN, not raise.
- Group names with `/`, spaces or `.` must not produce broken figure file names.
- PCoA with negative eigenvalues and with a genus that has missing codon rows: must drop or correct, and report how many.
- Pages deploy replacing the custom domain: the `CNAME` file must be in the deployed site.

---

## File Structure

| File | Responsibility |
|---|---|
| `bin/summary/make_figures.py` (modify) | Add `order` level, hue-shade colors, call new page functions, write figure index |
| `bin/summary/build_summary_tables.py` (modify) | New tables `gene_structure_genome`, `codon_genome`, `aa_genome` |
| `bin/summary/figs_structure.py` (create) | `page_gene_structure`: gene length, density, intron/exon, intergenic figures |
| `bin/summary/figs_composition.py` (create) | `page_composition`: codon and amino-acid figures |
| `bin/summary/stats.py` (create) | ANOVA/Kruskal tables and regression tables |
| `bin/summary/ordination.py` (create) | Genus representatives and PCoA + PERMANOVA runs; wraps `lib/ordination.py` |
| `bin/summary/figure_index.py` (create) | Build `figures_index.tsv` and a ranking table |
| `bin/summary/export_clade_json.py` (create) | Phase B: compact JSON for OJS |
| `docs/quarto_summary/*.qmd`, `_quarto.yml`, `_helpers.py` (modify/create) | New pages and the order tab |
| `.github/workflows/pages.yml`, `docs/quarto_summary/CNAME` (create) | Deployment |

---

### Task 1: Order level and phylum-hue shades

**Files:**
- Modify: `bin/summary/make_figures.py` (`LEVELS`, `assign_groups`, `PALETTE` use)
- Modify: `docs/quarto_summary/_helpers.py` (`LEVELS`)
- Test: `bin/summary/test_make_figures.py`

**Interfaces:**
- Produces: `LEVELS = ("phylum", "subphylum", "class", "order")`. `assign_groups(core, level)` still returns `(df, levels, colors)`. For `level="order"`, `grp = ORDER` with NULL → `"Unknown"`. `colors[g]` is a shade of the parent phylum's hue. New helper `shade(hex_color: str, i: int, n: int) -> str` (lightness varies, hue fixed).

- [ ] **Step 1: Write failing tests** in `test_make_figures.py`: `test_order_level_groups_exist` (`assign_groups(core, "order")` returns levels with no NULL), `test_shades_keep_parent_hue` (for two groups in one phylum: colors differ, hue within 0.02 in HSV via `colorsys`), `test_unknown_order_kept` (rows with NULL ORDER end up in `Unknown` or `Other`, total rows unchanged).
- [ ] **Step 2: Run** `python3.12 -m pytest -q bin/summary/test_make_figures.py -k "order or shade"`. Expected: FAIL.
- [ ] **Step 3: Implement** `shade` and the `order` branch of `assign_groups`. For each phylum, assign shade index by rank of group size within that phylum. Keep `OTHER` grey.
- [ ] **Step 4: Run** the full `test_make_figures.py`. Expected: PASS, including the existing tests (synthetic fixture uses `ORDER="o"`, so extend `synth` to vary `ORDER` and set some to NULL).
- [ ] **Step 5: Commit** `git commit -m "summary: order level and phylum-hue shades"`.

### Task 2: Per-genome structure and composition tables

**Files:**
- Modify: `bin/summary/build_summary_tables.py`
- Test: `bin/summary/test_build_summary_tables.py` (create; small DuckDB fixture built in the test)

**Interfaces:**
- Produces four parquet files written by `copy(con, name, sql, out)`:
  - `genome_all`: one row for EVERY genome in `v_genome_labels` (no `is_representative` filter). Columns: `species_prefix, PHYLUM, CLASS, ORDER, FAMILY, GENUS, SPECIES, is_representative, total_length_bp, gc_pct, busco_complete, busco_single, busco_duplicated, busco_fragmented, busco_missing (as available), N50_bp, contig_count, n_genes, mean_len, genes_per_mb, coding_fraction, median_intron_len_bp, median_intergenic_bp`. Read the BUSCO columns that exist in `busco_genome` and `asm_stats` and record which were found.
  - `gene_structure_genome`: `species_prefix, n_genes, median_gene_len_bp, median_cds_len_bp, median_protein_len_aa, genes_per_mb, coding_fraction, median_exons_per_gene, intron_frac_genes_with_intron, median_intron_len_bp, median_intergenic_bp`.
  - `codon_genome`: `species_prefix` plus 61 columns `codon_<AAA>`, frequencies summing to 1 within tolerance 1e-6.
  - `aa_genome`: `species_prefix` plus 20 columns `aa_<A>`.
- Consumes: DuckDB tables `gene_info`, `gene_proteins`, `gene_exons`, `gene_introns`, `gene_intergenic_distances`, `codon_freq`, `aa_freq`, and `reps` (existing temp table).

- [ ] **Step 1: Inspect** the real column names: `DESCRIBE` each source table through the dashboard DB (`db/BFD.duckdb`). Record them in the task notes. Gene tables are sorted by id: filter with `=` or a range, not `LIKE`. Group-by over the whole view can exceed 280 s: run in SLURM.
- [ ] **Step 2: Write failing test**: `genome_all` row count equals the `v_genome_labels` row count, including non-representatives; on a 3-genome fixture, `gene_structure_genome` has one row per genome, `genes_per_mb == n_genes / (total_length_bp/1e6)`, and `codon_genome` rows sum to 1.0 ± 1e-6.
- [ ] **Step 3: Run** it. Expected: FAIL.
- [ ] **Step 4: Implement** the three `copy(...)` calls. A genome with no introns gets `median_intron_len_bp = NULL`, not 0.
- [ ] **Step 5: Run** the test. Expected: PASS.
- [ ] **Step 6: Run on HPCC** with `sbatch` (job script in `workflows/summary_tables.sh`, `--partition epyc --cpus-per-task 16 --mem 64G`). Verify row counts: `genome_core` representatives vs `codon_genome` rows; report any difference. Output goes outside git.
- [ ] **Step 7: Commit.**

### Task 3: Structure and composition figures

**Files:**
- Create: `bin/summary/figs_structure.py`, `bin/summary/figs_composition.py`
- Modify: `bin/summary/make_figures.py` (`run_level` calls both, `Out` records topic)
- Test: `bin/summary/test_make_figures.py` (extend `synth` with the three new tables)

**Interfaces:**
- Produces: `page_gene_structure(core, levels, colors, structure, o)` and `page_composition(core, levels, colors, codon, aa, o)`. Both use `box_by_group`, `heat`, `stacked_by_group` from `make_figures`. `Out.save(name, topic=...)` gets a `topic` argument stored in `o.made` entries (name, topic, kind).
- Figure names: `S01_gene_len`, `S02_gene_density`, `S03_coding_fraction`, `S04_exons`, `S05_intron_len`, `S06_intron_frac`, `S07_intergenic`, `C01_codon_heatmap`, `C02_rscu_boxes`, `C03_aa_composition`, `C04_gc3`, `C05_aa_vs_gc`. Add more variants freely, because many figures are wanted. Keep names unique within a level and use only `[A-Za-z0-9_]`.

- [ ] **Step 1: Write failing test** `test_new_pages_make_figures`: run `mf.main(data, out, ("phylum", "order"))` on the extended fixture. Assert every listed name exists as PNG and SVG at both levels, and file names contain only allowed characters.
- [ ] **Step 2: Run.** Expected: FAIL.
- [ ] **Step 3: Implement** the two modules. One function per figure family, each wrapped so a failure in one figure logs and continues (`o.made` records skipped names).
- [ ] **Step 4: Run** all tests. Expected: PASS.
- [ ] **Step 5: Commit.**

### Task 4: Group tests, regressions and all-genome scatter plots

**Files:**
- Create: `bin/summary/stats.py`, `bin/summary/test_stats.py`
- Modify: `bin/summary/make_figures.py` (`run_level` writes `tables/<level>/T_anova.tsv`, `T_regression.tsv` and figures `R01_*`)

**Interfaces:**
- Produces: `group_tests(df: DataFrame, stats: list[str], group_col: str = "grp") -> DataFrame` with columns `stat, test, statistic, p_value, effect_size, n, n_groups`. Test is Kruskal-Wallis; effect size is epsilon² = H·(n+1)/(n²−1)… use `H/((n²−1)/(n+1))`. Add one-way ANOVA F and eta² as extra columns `anova_F, eta2`. Groups with fewer than 2 genomes are dropped.
- `pair_regressions(df, pairs: list[tuple[str, str]], log: bool = True) -> DataFrame` with columns `x, y, n, slope, intercept, r2, p_value`, via `scipy.stats.linregress` on log10 values when `log`.
- State on the page: regressions do not control for phylogeny.
- All-genome scatter family `A01`..`A12` (reserved; these use `genome_all`, not representatives and not genus reps): `A01_busco_vs_genome_size`, `A02_busco_vs_n_genes`, `A03_busco_vs_n50`, `A04_busco_vs_contigs`, `A05_busco_vs_gc`, `A06_busco_vs_gene_density`, `A07_genes_vs_genome_size`, `A08_intron_len_vs_gene_density`, `A09_gc_vs_genome_size`, `A10_busco_dup_vs_genome_size`, `A11_busco_missing_vs_contigs`, `A12_n50_vs_n_genes`. Signature: `scatter_all(df, x, y, colors_by_phylum, o, name, logx=False, logy=False)`. It draws a hexbin density layer plus points for phyla with fewer than 200 genomes in the phylum hue, and records Pearson r, Spearman r and n in the sidecar. The index marks them `scope=all_genomes`; other figures have `scope=representatives` or `scope=genus_reps`.

- [ ] **Step 1: Write failing tests**: `test_scatter_all_uses_every_genome` (fixture with 300 genomes of which 40 are representatives: the n recorded in the sidecar is 300); a synthetic frame with a known group shift gives `p_value < 1e-6` and `eta2 > 0.3`; a constant stat returns NaN effect size; a group of size 1 is dropped; `pair_regressions` on `y = 2x` returns `slope ≈ 2` and `r2 ≈ 1`.
- [ ] **Step 2: Run.** Expected: FAIL.
- [ ] **Step 3: Implement** both functions with `scipy.stats`.
- [ ] **Step 4: Run.** Expected: PASS.
- [ ] **Step 5: Wire** into `run_level`; make regression scatter figures per phylum-colored point (`R01_size_vs_genes`, `R02_density_vs_intron_len`, plus pairs from the topic stat lists).
- [ ] **Step 6: Commit.**

### Task 5: Genus representatives and PCoA

**Files:**
- Create: `bin/summary/ordination.py`, `bin/summary/test_ordination.py`
- Reuse: `lib/ordination.py` (read its PCoA/PERMANOVA signatures first and call them)

**Interfaces:**
- Produces: `pick_genus_reps(core: DataFrame) -> DataFrame`: one row per non-null `GENUS`, chosen by highest `busco_complete`, then highest `N50_bp`, then lowest `species_prefix`. Genomes with NULL genus are excluded and counted in the return `attrs["n_no_genus"]`.
- `run_ordination(reps, matrices: dict[str, DataFrame], out) -> dict`: for each of `codon`, `aa`, `structure` build a distance matrix (Bray-Curtis for codon and aa; Euclidean on standardized, log-transformed stats for structure), run PCoA (first 5 axes, share explained, negative eigenvalues reported), PERMANOVA by `PHYLUM` and by `ORDER` (999 permutations, seed 1). Write `tables/ordination/<name>_coords.tsv`, `<name>_permanova.tsv`, and figures `O01_<name>_pcoa_phylum`, `O02_<name>_pcoa_class` (shaded by phylum hue).

- [ ] **Step 1: Write failing tests**: `pick_genus_reps` returns unique genera and picks the BUSCO-best row on a fixture with ties; NULL genus rows are excluded and counted; PCoA on two well-separated synthetic clusters gives PERMANOVA `p < 0.01`; a genus with NaN codon rows is dropped and reported.
- [ ] **Step 2: Run.** Expected: FAIL.
- [ ] **Step 3: Implement.** Real size is about 1,671 genera: the matrix is small, so no HPCC job is needed.
- [ ] **Step 4: Run.** Expected: PASS.
- [ ] **Step 5: Run on real tables** and check the genus count against `select count(distinct GENUS) from species.parquet` (1,671 on 2026-11-17; the reps must not exceed it).
- [ ] **Step 6: Commit.**

### Task 6: Figure index, ranking and gallery data

**Files:**
- Create: `bin/summary/figure_index.py`, `bin/summary/test_figure_index.py`

**Interfaces:**
- Produces: `build_index(out_root: Path) -> DataFrame` scanning `figures/<level>/*.png`. It writes `figures_index.tsv` with `name, level, topic, kind, scope, n_genomes, n_groups, pooled_share, effect_size, score`. `effect_size` is joined from `T_anova.tsv` by statistic when the figure has one. `score = effect_size × (1 − pooled_share)`; figures with no effect size get an empty score and sort last.

- [ ] **Step 1: Write failing tests**: index has one row per PNG; `pooled_share` is in [0, 1]; a figure with no ANOVA row has empty `score`.
- [ ] **Step 2: Run.** Expected: FAIL.
- [ ] **Step 3: Implement.** Have `Out.save` write a small sidecar JSON next to each PNG with `topic, stat, n_genomes, n_groups, pooled_share` so the index does not parse names.
- [ ] **Step 4: Run.** Expected: PASS.
- [ ] **Step 5: Commit.**

### Task 7: Quarto pages and gallery

**Files:**
- Create: `docs/quarto_summary/07_gene_structure.qmd`, `08_composition.qmd`, `09_comparisons.qmd`, `10_ordination.qmd`, `11_gallery.qmd`
- Modify: `docs/quarto_summary/_quarto.yml` (navbar), `_helpers.py` (order tab, `gallery()` helper), `index.qmd` (overview)

**Interfaces:**
- Consumes: figure names from Tasks 3 to 5 and `figures_index.tsv`.
- Produces: `gallery(topic=None, level=None, top=None)` in `_helpers.py`, showing a filterable table (static HTML table with a client-side filter input) and thumbnails. Each content page starts with the phylum view, then tabs for subphylum, class, order.

- [ ] **Step 1: Write the pages** using `fig()` and `table()`. Keep captions factual; state n and the pooled groups. State the no-phylogeny caveat on `09_comparisons.qmd`. State the genus-representative rule on `10_ordination.qmd`.
- [ ] **Step 2: Render** locally: `module load quarto/1.7.31; QUARTO_PYTHON=<env>/bin/python quarto render docs/quarto_summary --output-dir $CLAUDE_JOB_DIR/tmp/site` (any scratch dir works). Expected: exits 0, all pages produce HTML.
- [ ] **Step 3: Check** that the render adds no untracked files to the repo (`git status --short`), and add `docs/quarto_summary/site_libs/`, `.quarto/`, `_site/` to `.gitignore` if it does.
- [ ] **Step 4: Check** the output size: `du -sh <site>` is below 600 MB.
- [ ] **Step 5: Commit.**

### Task 8: Deployment

**Files:**
- Create: `.github/workflows/pages.yml`, `docs/quarto_summary/CNAME`
- Modify: `docs/quarto_summary/_quarto.yml` (`resources: [CNAME]`)

**Interfaces:**
- Produces: workflow on push to `main` touching `docs/quarto_summary/**`: set up Python 3.12, `pip install pandas matplotlib nbformat nbclient ipykernel pyyaml tabulate`, set up Quarto 1.7.x via `quarto-dev/quarto-actions/setup`, `quarto render docs/quarto_summary`, upload `_site` with `actions/upload-pages-artifact`, deploy with `actions/deploy-pages`. Permissions: `pages: write`, `id-token: write`. `CNAME` content: `bfd.fungalgenomes.org`.
- Figures and small tables must be committed under `docs/quarto_summary/figures` and `tables` (current state: already tracked; verify with `git ls-files | head`).

- [ ] **Step 1: Test locally**: after render, `test -f <site>/CNAME && cat <site>/CNAME` prints the domain.
- [ ] **Step 2: Validate** the workflow YAML (`python3.12 -c "import yaml,sys;yaml.safe_load(open('.github/workflows/pages.yml'))"`).
- [ ] **Step 3: Commit and open a PR.** Do not merge or push to `main` without the user's approval.
- [ ] **Step 4: User steps (not automatable here):** in repo Settings → Pages set Source to GitHub Actions, set the custom domain, enable HTTPS. Report the check once the DNS resolves.

### Task 9: Phase B start: JSON export and one OJS page

**Files:**
- Create: `bin/summary/export_clade_json.py`, `bin/summary/test_export_clade_json.py`, `docs/quarto_summary/12_explore.qmd`

**Interfaces:**
- Produces: `export_clade_json(core: DataFrame, structure: DataFrame, out: Path) -> Path` writing `docs/quarto_summary/data/genomes.json` (and `.json.gz` is not served; keep plain JSON under 5 MB). Records: `id, phylum, subphylum, class, order, genus, genome_mb, n_genes, gc, genes_per_mb, median_intron_len_bp, median_intergenic_bp, coding_fraction`. Round floats to 3 digits. Null taxonomy → `"Unknown"`.
- `12_explore.qmd`: an `{ojs}` page that loads `data/genomes.json`, with selectors for x, y and color level (phylum, class, order) and a phylum filter. Colors follow the Task 1 shade scheme: export a `colors` map into the same JSON.

- [ ] **Step 1: Write failing tests**: output parses as JSON; record count equals the number of rows; file is under 5 MB for the real table; no NaN literal appears in the file.
- [ ] **Step 2: Run.** Expected: FAIL.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run.** Expected: PASS.
- [ ] **Step 5: Render** the page and open `12_explore.html` check: the JSON loads (`resources: data/genomes.json` in `_quarto.yml` so it deploys).
- [ ] **Step 6: Commit.**

---

## Not in this plan

- Phase C (DuckDB-WASM with an LLM). It cannot run on GitHub Pages as a data host. Needs its own design and a hardware study. See `docs/WASM_PLAN.md`.
- Fixing stale tables listed in `WASM_PLAN.md` section 3.

## Self-review notes

- Spec coverage: all-genome scatter plots including BUSCO (T2 `genome_all`, T4 `A01`..`A12`), levels (T1), topics and tables (T2, T3), ANOVA/regression (T4), PCoA with genus reps (T5), figure index, ranking and gallery (T6, T7), hosting (T8), phase B (T9). Phase C is excluded by decision.
- Names used across tasks: `assign_groups`, `Out.save(topic=)`, `group_tests`, `pair_regressions`, `pick_genus_reps`, `run_ordination`, `build_index`, `export_clade_json`: each defined once.
- Open: column names in Task 2 depend on a schema check (Step 1 of Task 2). Task 2 is the only task that needs the real database.
