# BFD TF Survey Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Enumerate transcription factors across all BFD genomes, test the lineage claims of Lax et al. 2025 at that scale, track WC-2, Zap1-like and bZip5 across fungi, and search for TF-family expansions and contractions the paper did not cover. Results appear as four new Quarto pages.

**Architecture:** Extend the existing summary pipeline (`bin/summary/`). A SLURM table builder writes small per-genome TF tables. New figure/statistics modules plug into `make_figures.main` like the existing ones (tables optional, suppress list applied at read time, outputs through `Out`). Ortholog tracking is a separate small SLURM workflow in `bin/tf/` that writes tables and figures into the same output tree. Outputs go into the release asset, not git.

**Tech Stack:** Python 3.12 (pandas, numpy, scipy, matplotlib, pyarrow, duckdb, biopython 1.88, pytest); HPCC modules `MMseqs2/17-b804f`, `mafft/7.505`, `fasttree/2.1.11`, `hmmer/3.4`, `iqtree/3.0.1` (final trees only); Quarto 1.7.31; `gh` for the release.

**Spec:** `docs/superpowers/specs/2026-11-17-bfd-tf-survey-design.md`. Paper: `resources/s41467-025-65177-1.pdf`.

## Global Constraints

- Branch `tf-survey` is stacked on `quarto-site` (PR #18). No push, PR or release without the user's approval.
- Test command: `/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/atlas_env/venv/bin/python -m pytest -q -W error::UserWarning bin/summary/ bin/tf/` (the existing suite takes about 13 minutes; run focused files while iterating and the full suite once before each commit series).
- Quarto 1.7.31 renders 0 files with exit 0 under a dot-folder path (`.claude/...`): always render from a `/tmp` copy.
- SLURM: small requests (`-c 2 --mem 24G` for figure/stat jobs; `-c 8 --mem 32G` for the mmseqs job); `$SCRATCH:?` for temp files only; final outputs on `/bigdata`; NEVER `BASH_SOURCE`; write the SLURM log with `-o %x.%j.log` (submit directory); jobs may queue for hours (account CPU limit): do not block on them for more than a few minutes, report job ids.
- Genome key is `species_prefix` (= LOCUSTAG). Never derive it with `left(protein_id, 8)` (three LOCUSTAGs have 9 characters): map protein to genome through `gene_proteins`.
- The suppress list (`config/summary_suppress.tsv`, module `bin/summary/suppress.py`) is applied at read time to every new table that has `species_prefix`.
- All tests use group values at genus level (one value per genus, from genus representatives) and Benjamini-Hochberg correction. Phylogeny is not otherwise controlled; every page says so. No regulation or binding claim for any genome except as quoted from the paper.
- One global phylum colour map (`phylum_color_map` in `bin/summary/make_figures.py`) is used by every figure.
- Compress large text output (zstd or gzip); Python readers accept compressed input where they read large text.
- Plain, short-sentence technical English in pages and docs.

## Review Focus

- A Pfam accession in the family map that is missing from `pfam_names`: the builder must fail loudly, not drop the family.
- A genome with zero TF proteins, a zero-size proteome, or one family only: fractions are NaN (not 0, not a crash); Shannon index is 0 for one family and NaN for none.
- A family present in fewer than 5 genera, a lineage with fewer than 5 genera, a genus with all-NaN values: tests are skipped and reported, never raise.
- The bootstrap and permutation steps are deterministic (fixed seed) and the same on every run.
- A seed protein that cannot be fetched (network error, not found): clear error, cache used on rerun, run continues with the other seeds and reports which are missing.
- Orthologs: a genus representative without the target domain is reported as absent, not dropped; domain-architecture rules are applied after the search, not before.
- Suppressed genomes never appear in any TF table or figure.

---

## File Structure

| File | Responsibility |
|---|---|
| `config/tf_family_map.tsv` (create) | Pfam accession → BFD family → paper family, with notes |
| `config/tf_claims.tsv` (create) | One row per paper claim: id, text, family, metric, group A, group B, direction |
| `config/tf_seeds.tsv` (create) | Seed proteins for ortholog tracking: name, source, key, target domain, required architecture |
| `bin/summary/build_tf_tables.py` (create) | `tf_family_genome`, `tf_summary_genome` for ALL genomes (DuckDB over `tf_inventory` + `pfam`) |
| `bin/summary/tf_data.py` (create) | Loaders and genus aggregation shared by the TF modules |
| `bin/summary/figs_tf.py` (create) | Landscape figures `T01..T08` by level (representatives) and the all-genome scatter |
| `bin/summary/tf_claims.py` (create) | Claims table and forest figure |
| `bin/summary/tf_expansion.py` (create) | Scaling model, lineage tests, presence/absence flags, heatmap |
| `bin/tf/seeds.py`, `bin/tf/ortho_search.py`, `bin/tf/ortho_tree.py`, `bin/tf/ortho_report.py` (create) | Seed fetch, search, tree, tables and figures |
| `workflows/tf_tables.sh`, `workflows/tf_orthologs.sh` (create) | SLURM scripts |
| `bin/summary/make_figures.py` (modify) | Call the TF modules when the TF tables exist |
| `docs/quarto_summary/13_tf_landscape.qmd` … `16_tf_expansions.qmd`, `_quarto.yml`, `_helpers.py`, `06_data_notes.qmd` (create/modify) | Pages |

---

### Task 1: TF family map and per-genome TF tables

**Files:**
- Create: `config/tf_family_map.tsv`, `bin/summary/build_tf_tables.py`, `bin/summary/test_build_tf_tables.py`, `workflows/tf_tables.sh`

**Interfaces:**
- Produces: `load_family_map(path: Path) -> DataFrame` with columns `pfam_acc, family, paper_family, note` (tab-separated file; `#` comments and blank lines ignored; duplicate accessions are an error).
- Produces: `validate_accessions(fmap: DataFrame, pfam_names: DataFrame) -> None`, raising `ValueError` that lists every accession missing from `pfam_names` (`pfam_acc` column, versions stripped).
- Produces: `build_tf_tables(con, tables_dir: Path, summary_dir: Path, fmap: DataFrame, out: Path) -> dict` writing `out/tf_family_genome.parquet` (`species_prefix, family, n_proteins`, one row per genome and family with n ≥ 1, ZSTD) and `out/tf_summary_genome.parquet` (`species_prefix, n_proteins_total, n_tf_proteins, tf_fraction, n_families, shannon`), for EVERY genome in `genome_all` (no representative filter; genomes with no TF rows still get a summary row with `n_tf_proteins = 0`).
- Family set: the 21 existing BFD families from `nextflow/assets/fungal_tf_pfam_domains.csv` plus E2F_TDP (PF02319), Tub (PF01167), RFX_DNA_binding (PF02257), WRKY (PF03106), CP2 (PF04516), CBFD_NFYB_HMF (PF00808), Copper-fist (PF00649), TEA (PF01285), KilA-N (PF04383); STAT: look up the STAT Pfam accessions in `pfam_names` by name (`STAT_bind`, `STAT_alpha`, `SH2`-containing STAT entries); add them if they exist, otherwise list STAT as `unmapped` in the notes column of a row with an empty `pfam_acc` (the loader keeps such rows out of the counts but the data-notes page shows them).
- Counting rule: proteins with at least one domain of the family; existing 22 Pfams come from `tf_inventory` (it has `species_prefix`, `protein_id`, `pfam_acc`, `tf_family`), the added ones from `pfam` filtered by accession and mapped to `species_prefix` through `gene_proteins(species_prefix, protein_id)`. Read `nextflow/bin/merge_tf_domains.py` and apply exactly the same domain gate (e-value/score) to the added Pfams; state the gate in the module docstring.
- Proteome size: `n_proteins_total` = distinct proteins in `gene_proteins` per genome. Shannon: natural log of the family distribution (proteins per family / `n_tf_proteins`), 0 for one family, NaN for none.

- [ ] **Step 1: Write the failing tests** (`test_build_tf_tables.py`, tiny parquet fixtures): `test_family_map_rejects_duplicate_accessions`, `test_validate_accessions_lists_all_missing` (two missing accessions named in the message), `test_tables_cover_every_genome_including_no_tf` (3-genome fixture: genome with no TF rows has `n_tf_proteins == 0` and `tf_fraction == 0.0` when the proteome is non-empty, NaN when the proteome is empty), `test_counts_proteins_not_domains` (a protein with two Zn_clus domains counts once), `test_added_pfam_mapped_through_gene_proteins` (a protein id whose genome is a 9-character prefix maps correctly), `test_shannon_values` (one family → 0.0; two equal families → ln 2 ≈ 0.693147), `test_same_gate_as_pipeline` (a hit below the gate is excluded).
- [ ] **Step 2: Run** `pytest -q bin/summary/test_build_tf_tables.py`. Expected: FAIL (module missing).
- [ ] **Step 3: Implement** the three functions and a CLI `python bin/summary/build_tf_tables.py <tables_dir> <summary_tables_dir> <out_dir>` (paths resolved from `Path(__file__)`; `--map` default `config/tf_family_map.tsv`; ZSTD parquet). Create `config/tf_family_map.tsv` with the families above, `paper_family` set from the paper's wording (for example `Zn cluster`, `C2H2/CCHC/CH/C5HC2`, `HTH/Homeodomain-like`, `HLH`, `GATA Zn finger`, `bZip`, `HMG-box`, `E2F`, `Tubby`, `RFX`, `WRKY`, `HSF`, `MADS`, `Forkhead`, `CP2`, `CBF`, `Copper fist`, `APSES`, `TEA`) and a `note` explaining each non-obvious mapping.
- [ ] **Step 4: Run** the tests. Expected: PASS.
- [ ] **Step 5: Write `workflows/tf_tables.sh`** (SLURM `-c 8 --mem 32G --time 2:00:00 -J tf_tables -o %x.%j.log`, `ROOT` from `${ROOT:-${SLURM_SUBMIT_DIR:?}}`, DuckDB temp in `$SCRATCH`, output to `/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/summary_tables`, `bash -n` clean, mode 755). Submit it from the worktree root and record the job id. Expected on completion: both parquet files exist; row count of `tf_summary_genome` equals the `genome_all` row count; `n_tf_proteins` is 0 for the genomes with no proteome; every accession in the map was found; STAT handling reported.
- [ ] **Step 6: Commit.**

### Task 2: Landscape figures by lineage

**Files:**
- Create: `bin/summary/tf_data.py`, `bin/summary/figs_tf.py`, `bin/summary/test_figs_tf.py`
- Modify: `bin/summary/make_figures.py` (`run_level`, `main`)

**Interfaces:**
- Consumes: `Out.save(name, topic=...)`, `box_by_group`, `heat`, `stacked_by_group`, `assign_groups`, `phylum_color_map`, `figs_common.safe`, `suppress.drop_suppressed` (read-time filter in `main`'s `rd`).
- Produces: `page_tf_landscape(core, levels, colors, fam, summ, o) -> None` where `fam` = `tf_family_genome`, `summ` = `tf_summary_genome`, `core` = the level frame from `assign_groups` (representatives). Figures (PNG+SVG, topic `tf`): `T01_tf_fraction` (box of `tf_fraction` by group), `T02_tf_composition` (stacked mean fraction of TF proteins by family), `T03_tf_family_heatmap` (group × family, z-score of the per-genome count per Mb of proteome... use per-genome `n_proteins / n_proteins_total`, mean by group, z-scored by family), `T04_tf_rare_presence` (heatmap, fraction of genomes with ≥ 1 protein, for E2F, Tubby, RFX, WRKY, CP2, STAT if present), `T05_tf_shannon` (box), `T06_tf_vs_proteome` (all-genome scatter via `figs_scatter.scatter_all`, x = proteome size, y = `n_tf_proteins`, both log; scope `all_genomes`).
- Produces: `tf_data.load_tf(data_dir: Path, suppress) -> tuple[DataFrame, DataFrame]`, `tf_data.genus_table(core_all, fam, summ, busco_min=None) -> DataFrame` (one row per genus: the best-BUSCO genome with complete data, columns `genus, species_prefix, PHYLUM, CLASS, ORDER, busco_complete, n_proteins_total, n_tf_proteins, tf_fraction, shannon` plus one `fam_<family>` column per family with the protein count, 0 where absent; excludes genomes with `busco_complete < busco_min` when given; NULL genus excluded and counted in `attrs["n_no_genus"]`).
- `make_figures.main`: when `tf_family_genome.parquet` and `tf_summary_genome.parquet` exist, `run_level` calls `page_tf_landscape`; when absent, print one note and skip (existing runs keep working). Level pages use representatives like the other pages.

- [ ] **Step 1: Write the failing tests** (fixture frame, 5 phyla, 3 levels): `test_all_tf_figures_written` (names T01–T06 exist as PNG+SVG at phylum and order level), `test_main_skips_without_tf_tables` (existing fixture without TF tables: `main` runs, prints the skip note), `test_genus_table_one_row_per_genus_and_best_busco` (tie rules as `pick_genus_reps`: busco, then N50, then lowest id), `test_genus_table_busco_filter` (a genus whose only genome is below the cutoff disappears and is counted), `test_nan_and_empty_groups_do_not_crash` (a group with one genome, a family absent everywhere), `test_suppressed_genome_absent_from_tf_outputs`, `test_zero_tf_genome_is_nan_fraction_safe`.
- [ ] **Step 2: Run** `pytest -q bin/summary/test_figs_tf.py`. Expected: FAIL.
- [ ] **Step 3: Implement** `tf_data.py` and `figs_tf.py`; each figure inside `safe()`; names fixed constants; titles state the unit ("proteins with a TF domain / all proteins"). Wire `main`/`run_level`.
- [ ] **Step 4: Run** the focused tests, then `pytest -q bin/summary/test_make_figures.py bin/summary/test_figure_index.py` (the index must pick up the new figures automatically; add `T`-figure stats to `figure_index` only if the sidecar path needs it).
- [ ] **Step 5: Commit.**

### Task 3: The paper's claims at BFD scale

**Files:**
- Create: `config/tf_claims.tsv`, `bin/summary/tf_claims.py`, `bin/summary/test_tf_claims.py`

**Interfaces:**
- Consumes: `tf_data.genus_table`.
- Produces: `config/tf_claims.tsv` columns `claim_id, text, family, metric, level, group_a, group_b, paper_direction` where `metric` is `fraction` (family proteins / all proteins), `count` or `prevalence` (fraction of genera with ≥ 1), `group_a` and `group_b` are phylum or lineage names or the sets `Dikarya` (Ascomycota + Basidiomycota) and `EDF` (all other phyla except Microsporidia; the file's header comment says this is the controller's reading of the paper), and `paper_direction` is `a>b`, `a<b` or `absent_in_a`. One row per claim listed in the spec section "What the paper claims that we can test".
- Produces: `evaluate_claims(genus_df: DataFrame, claims: DataFrame, n_boot: int = 2000, seed: int = 1) -> DataFrame` with columns `claim_id, text, family, metric, group_a, group_b, n_genera_a, n_genera_b, value_a, value_b, effect, ci_lo, ci_hi, p, p_adj, observed_direction, paper_direction, agrees`. `value_*` = median per-genus metric (or prevalence); `effect` = rank-biserial correlation (Mann-Whitney U based) for `fraction`/`count`, and difference in prevalence for `prevalence`; CI = percentile bootstrap over genera (`n_boot` resamples within each group, seeded); `p` = Mann-Whitney (Fisher exact for prevalence); `p_adj` = Benjamini-Hochberg over all claims; `agrees` = `yes` when the sign of the effect matches `paper_direction` and the CI excludes 0, `mixed` when the sign matches but the CI includes 0, `no` when the sign is opposite. A group with fewer than 5 genera gives NaN values and `agrees = "not tested"` with the reason in a `note` column.
- Produces: `page_tf_claims(genus_df, claims_table, o)` writing `tables/tf/T10_tf_claims.tsv` and the forest figure `T10_tf_claims_forest` (effect with CI per claim, coloured by `agrees`).

- [ ] **Step 1: Write the failing tests**: `test_rank_biserial_known_answer` (groups [1,2,3] vs [4,5,6]: effect −1.0; [1,2,3] vs [3,4,5] with one tie: hand-computed), `test_bootstrap_is_deterministic` (same seed, same CI; different seed, different CI), `test_bh_adjustment_known_answer` (p = [0.01, 0.04, 0.03] → [0.03, 0.04, 0.04]), `test_agrees_logic` (sign + CI cases: yes, mixed, no), `test_group_with_fewer_than_5_genera_not_tested`, `test_prevalence_claim_uses_fisher` (counts hand-checked), `test_edf_and_dikarya_sets_resolved`, `test_missing_family_column_reported_not_raised`.
- [ ] **Step 2: Run.** Expected: FAIL.
- [ ] **Step 3: Implement** with `scipy.stats` and `numpy.random.default_rng(seed)`. Write the claims file.
- [ ] **Step 4: Run** the tests. Expected: PASS.
- [ ] **Step 5: Wire** `page_tf_claims` into `make_figures.main` (after the landscape page; BUSCO cutoff 80 and a sensitivity table without the cutoff, `T10_tf_claims_nofilter.tsv`).
- [ ] **Step 6: Commit.**

### Task 4: Expansions and contractions

**Files:**
- Create: `bin/summary/tf_expansion.py`, `bin/summary/test_tf_expansion.py`

**Interfaces:**
- Consumes: `tf_data.genus_table`.
- Produces: `fit_scaling(genus_df: DataFrame, family: str) -> DataFrame` adding `resid` = residual of `log1p(count)` regressed on `log(n_proteins_total)` (ordinary least squares over all genera with a proteome size; NaN where either is missing), plus slope/intercept in `attrs`.
- Produces: `lineage_tests(genus_df, families: list[str], level: str = "ORDER", min_genera: int = 5, min_agree: int = 3) -> DataFrame` columns `family, lineage, parent, n_genera, n_rest, median_resid, effect, p, p_adj, direction, n_agree, flagged`; each lineage is compared with the other genera of the same phylum (Mann-Whitney on residuals; effect = rank-biserial; BH over all family-by-lineage tests; `n_agree` = genera in the lineage whose residual has the same sign as the median; `flagged` when `p_adj < 0.05`, `|effect| >= 0.5` and `n_agree >= min_agree`).
- Produces: `presence_flags(genus_df, families, level="ORDER", min_genera=5) -> DataFrame` columns `family, lineage, parent, global_prevalence, lineage_prevalence, flag` where `flag` is `lost` when global prevalence ≥ 0.8 and lineage prevalence ≤ 0.2, `gained` for the reverse (global ≤ 0.2, lineage ≥ 0.8), else empty.
- Produces: `page_tf_expansion(genus_df_busco, genus_df_all, o)` writing `tables/tf/T11_tf_expansion_hits.tsv`, `T11_tf_presence_flags.tsv`, `T11_tf_expansion_nofilter.tsv` (sensitivity: no BUSCO filter), figures `T11_tf_effect_heatmap` (family × lineage effect, flagged cells outlined, rows grouped by phylum colour), `T12_tf_hit_boxes` (box of residual for up to 12 top flagged hits against the rest of the phylum).

- [ ] **Step 1: Write the failing tests**: `test_fit_scaling_recovers_known_slope` (synthetic counts = 2 × proteome^0.8, slope ≈ 0.8, residuals ≈ 0), `test_lineage_test_flags_planted_expansion` (one order with +1.5 SD residuals over 8 genera inside a phylum of 40 → flagged), `test_no_flag_below_min_genera` (4 genera → not tested, note), `test_bh_across_all_tests` (two families × three lineages; adjusted values monotone), `test_presence_flags_lost_and_gained`, `test_constant_family_gives_nan_not_error` (a family with count 0 everywhere), `test_deterministic_output_order` (sorted by `p_adj`, then family, then lineage), `test_sensitivity_run_differs_when_low_busco_genera_planted`.
- [ ] **Step 2: Run.** Expected: FAIL.
- [ ] **Step 3: Implement** with numpy/scipy; compute OLS with `numpy.polyfit`; tests reuse BH helper from Task 3 (`tf_claims.bh_adjust`, move it into `bin/summary/stats.py` as `bh_adjust(p) -> ndarray` and import in both; test it there).
- [ ] **Step 4: Run** the tests. Expected: PASS.
- [ ] **Step 5: Wire** into `make_figures.main` (guarded like the others; hits are labelled "hypothesis" in titles).
- [ ] **Step 6: Commit.**

### Task 5: Ortholog tracking (WC-2, Zap1-like, bZip5)

**Files:**
- Create: `config/tf_seeds.tsv`, `bin/tf/__init__.py`, `bin/tf/seeds.py`, `bin/tf/ortho_search.py`, `bin/tf/ortho_tree.py`, `bin/tf/ortho_report.py`, `bin/tf/test_tf_ortho.py`, `workflows/tf_orthologs.sh`

**Interfaces:**
- `config/tf_seeds.tsv` columns `target, seed_name, source, key, taxon_id, target_pfam, required_pfam, note`. `source` is `uniprot` (key = gene symbol, `taxon_id` the organism: *S. cerevisiae* 559292 for `ZAP1` and `CAD1`, *N. crassa* 367110 for `wc-2` locus NCU00902) or `bfd` (key = `species_prefix` plus a protein selector `pfam:GATA+PAS` meaning "all proteins in that genome with both domains"; used for the BFD *R. microsporus* ATCC 11559 `FBF4BF29`, *M. lusitanicus* `F750FE6D` and *P. blakesleeanus* `F43EADC6` white-collar candidates). Targets: `WC-2` (target_pfam PF00320 GATA; required PF00989 PAS_*; accept any PAS Pfam: list the accessions in `required_pfam` separated by `|`), `Zap1-like` (PF00096 C2H2 with ≥ 4 domains), `bZip5` (PF00170 bZIP_1). Verify every accession against `pfam_names`.
- `seeds.fetch_seeds(seeds: DataFrame, cache_dir: Path, fetch=urllib/requests) -> DataFrame` (columns `target, seed_name, accession, sequence, status`): UniProt REST (`https://rest.uniprot.org/uniprotkb/search?query=gene_exact:<key>+AND+organism_id:<taxon>&format=fasta&size=1`), cached by key; network or not-found errors give `status = "missing: <reason>"`, never an exception; `fetch` is injectable so tests use a stub. A `bfd` seed is resolved from `tf_inventory`/`pfam` and `seqdb/proteins`.
- `ortho_search.select_candidates(genus_reps: DataFrame, pfam_hits: DataFrame, target: str, seeds: DataFrame) -> DataFrame` (all proteins of genus representatives carrying the target Pfam, with their domain architecture string); `ortho_search.search(seed_fasta: Path, cand_fasta: Path, out: Path, scratch: Path, min_pident=30.0, min_cov=0.5) -> DataFrame` (MMseqs2 `easy-search`; thresholds written into the output table header/attrs); `ortho_search.architecture_ok(arch: str, required: list[str]) -> bool`.
- `ortho_tree.build_tree(fasta: Path, out_prefix: Path, scratch: Path) -> Path` (`mafft --auto`, trim columns with > 50% gaps, `FastTree -lg`; returns the Newick path); `ortho_tree.draw_tree(newick: Path, meta: DataFrame, colors: dict, out_png: Path)` (Bio.Phylo + matplotlib, tip points coloured by the global phylum colour map, anchor species marked).
- `ortho_report.summarise(hits: DataFrame, genus_reps: DataFrame) -> DataFrame` columns `target, lineage, parent, n_genera, n_with_ortholog, prevalence, median_copies, architecture_counts` for phylum and order; written to `tables/tf/T20_ortho_<target>.tsv`; figure `T20_ortho_<target>_presence` (order × presence heatmap) and `T20_ortho_<target>_tree`.
- `workflows/tf_orthologs.sh`: SLURM (`-c 8 --mem 32G --time 3:00:00 -J tf_ortho -o %x.%j.log`), loads modules in separate lines (not piped), runs seeds → candidates → search → tree → report per target, writes to `$OUT` (default the staging figures dir), temp in `$SCRATCH`.

- [ ] **Step 1: Write the failing tests** (`test_tf_ortho.py`, no network, no external binaries except where skipped with a clear reason): `test_seed_fetch_uses_cache_and_reports_missing` (stub fetch: success cached; second call does not call fetch; a failing key gives `status` starting with `missing`), `test_architecture_ok_rules` (GATA+PAS present → true; GATA only → false; `>=4 C2H2` rule), `test_select_candidates_includes_genera_without_hits_as_absent` (genus rep without the domain is kept in the report as absent), `test_report_prevalence_and_copies_known_answer`, `test_thresholds_recorded_in_output`, `test_build_tree_skips_with_clear_message_when_mafft_missing` (monkeypatch `shutil.which`), `test_search_skipped_when_mmseqs_missing`, and an integration test marked to skip unless `mmseqs`, `mafft` and `FastTree` are on PATH that runs the chain on 12 toy sequences and gets a Newick with 12 tips.
- [ ] **Step 2: Run.** Expected: FAIL.
- [ ] **Step 3: Implement** the modules and the SLURM script (`bash -n` clean). Pull genus representatives from `genome_all` (best BUSCO per genus, as in the ordination) and proteins from `seqdb/proteins` through `samtools faidx` or the existing helper (look at `bin/bfd_scratch.py`/`lib/` first and reuse).
- [ ] **Step 4: Run** the focused tests. Expected: PASS (integration test passes with modules loaded, otherwise skips with its reason).
- [ ] **Step 5: Submit** `workflows/tf_orthologs.sh` for the real data (after Task 1's tables exist); record the job id; do not wait more than a few minutes; when it finishes check: seeds resolved (status table), candidate counts per target, thresholds, tree tip counts, presence table for `WC-2` in Mucoromycota, Ascomycota, Basidiomycota, Chytridiomycota.
- [ ] **Step 6: Commit.**

### Task 6: Pages

**Files:**
- Create: `docs/quarto_summary/13_tf_landscape.qmd`, `14_tf_contrasts.qmd`, `15_tf_orthologs.qmd`, `16_tf_expansions.qmd`
- Modify: `docs/quarto_summary/_quarto.yml`, `_helpers.py`, `06_data_notes.qmd`, `bin/summary/test_site_pages.py`

**Interfaces:**
- Consumes: figure names `T01..T06`, `T10_*`, `T11_*`, `T12_*`, `T20_*`; tables `tables/tf/T10_tf_claims.tsv`, `T11_*`, `T20_ortho_*`; existing helpers `fig`, `single`, `table`, `flat_table`, `read_flat`, `gallery`.
- Required statements: `14` says: effects are per genus (best-BUSCO genome, BUSCO ≥ 80%), phylogeny is not otherwise controlled, p-values are Benjamini-Hochberg adjusted, and "agrees with the paper" only compares the sign; it names the paper (Lax et al. 2025, Nature Communications 16:10277) and that the paper used InterPro-based families and 62 species. `15` says: seed proteins and sources, MMseqs2 thresholds, the architecture rule, that similarity and tree position are evidence for orthology, not proof, and that function is not inferred; anchor species are marked. `16` says: hits are hypotheses; annotation completeness affects counts; the sensitivity table is shown; effect cutoffs (`p_adj < 0.05`, `|effect| ≥ 0.5`, ≥ 3 genera agreeing). `06_data_notes` lists the TF family map including unmapped paper families (read from `config/tf_family_map.tsv` copy under `data/` with a sync test, like the suppress list) and the build commands for the TF tables.
- Pages degrade gracefully when a figure or table has not been built ("has not been built yet" notes); the render of the committed tree must exit 0 with the asset missing for these pages.

- [ ] **Step 1: Write the failing tests** (`test_site_pages.py`): navbar hrefs exist for the four pages; each required statement appears in its page by substring; the family-map copy equals the config file; the pages render with no TF tables (graceful notes); `gallery()` still lists TF figures when present in `figures_index.tsv`.
- [ ] **Step 2: Run.** Expected: FAIL.
- [ ] **Step 3: Write** the pages (figures through `fig()`/`single()`, tables through `flat_table`, captions factual with n read from sidecars) and update navbar after "Explore".
- [ ] **Step 4: Run** the tests, then render from a `/tmp` copy (fixture figures for the preview, as in the earlier site work): exit 0, 17 HTML pages.
- [ ] **Step 5: Commit.**

### Task 7: Real run, new release asset, PR

**Files:**
- Modify: `docs/quarto_summary/figures_release.json` (by the release script)

- [ ] **Step 1:** Confirm the SLURM jobs from Tasks 1 and 5 finished (or run Task 2–4 figure code on the login node at low priority as before if the queue is stuck: `nice -n 10`, one thread, output to the staging dir). Run `python bin/summary/make_figures.py <summary_tables> <staging_out>`; expected exit 0, the index lists the TF figures, 0 skipped except documented ones.
- [ ] **Step 2:** Look at one figure per family of figures (T01, T03, T04, T10 forest, T11 heatmap, a tree) and the claims table; list anything unreadable or surprising (for example a claim with `agrees = no`: report it, do not hide it).
- [ ] **Step 3 (user-gated, outward-facing):** `workflows/release_site_assets.sh --date <new date> --upload` creates a NEW release (new dated tag; never `--clobber` the published one); commit the new manifest; ask the user before running it.
- [ ] **Step 4:** Update the stacked PR text; no push or merge without the user's approval.

---

## Self-review notes

- Spec coverage: A (T1, T2, T3), B (T5), C (T4), pages (T6), release and real run (T7). Anti-goals: no verdict column beyond the sign (T3 `agrees`), BUSCO filter + sensitivity (T2–T4).
- Names used across tasks: `genus_table`, `load_tf`, `evaluate_claims`, `bh_adjust` (moved into `bin/summary/stats.py` in Task 4 and imported by Task 3 only after Task 4 moves it — until then Task 3 has its own `bh_adjust` in `tf_claims.py`; Task 4 replaces it with the shared one and keeps the tests green).
- Open: STAT Pfam accessions (Task 1 looks them up); the paper's supplementary tables (optional input for Task 5; seeds table can add rows with `source = paper` later).
