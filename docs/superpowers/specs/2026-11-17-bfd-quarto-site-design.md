# BFD Quarto site at bfd.fungalgenomes.org: design

Date: 2026-11-17. Status: DRAFT, awaiting approval. Repo: `stajichlab/BFD_dashboard`.

## Intent

Publish a browsable site of genome, gene and composition statistics for the BFD fungal genomes (23,680 rows in `species`). Readers see a high-level view first. They then drill into phylum, order and lower groups. The user wants many figures first and a strategy to choose the best ones afterwards.

## Goals

- Topics: genome stats, gene stats, gene length, codon usage, amino-acid usage, gene density, gene structure (introns, exons, intergenic distance).
- Comparisons between topics: ANOVA and regression, and PCoA with related ordinations.
- Levels: phylum, subphylum, class, order (new). Subgroups use shades of the parent phylum's hue.
- Generate a large figure set with one naming scheme, plus a figure index.
- A gallery page shows all figures. A selection strategy promotes the best to main pages.
- Host on GitHub Pages at `bfd.fungalgenomes.org` (DNS CNAME to `stajichlab.github.io` already exists).

## Non-goals (this phase)

- No live filtering or browser-side queries. That is phase B (Observable JS) and phase C (DuckDB-WASM).
- No new annotation runs. Use existing tables only.
- No ranking by an algorithm alone. The user makes the final choice.

## Anti-goals

- Limiting scatter plots to one genome per genus. A reserved family of scatter plots (BUSCO against size, genes, N50, contigs, GC, gene density) uses ALL genomes. PCoA stays on genus representatives.
- Publishing numbers from stale or incomplete tables without a caveat. `targetP` (125 genomes) and `infernal_rfam` (14) are excluded. Tables built before the 2026-11-16 fixes are re-checked first.
- Committing large data to the Pages repo history. Pages limit is 1 GB per site. Current figures are 132 MB.

## Constraints

- Quarto `module load quarto/1.7.31`. It needs a Python env with `nbformat`, `nbclient`, `ipykernel`, `pandas`, `tabulate`, `matplotlib`. A test render of `docs/quarto_summary` succeeded in 38 s with such an env (not yet recorded in the repo).
- Heavy computation runs on HPCC through SLURM, on genus-level or smaller tables. Quarto only displays results.
- The pipeline code must not touch files under a live Nextflow run.
- Genomes lacking taxonomy: 341 have no genus and 300 have no order. They are kept at the levels where taxonomy exists and labelled "unassigned" elsewhere.

## Approaches considered

- **A. Static pipeline extension (chosen by the user).** Extend `bin/summary/make_figures.py` and `taxonomy.py`. Quarto shows pre-built PNG/SVG and tables.
- **B. Observable JS pages on per-clade JSON (later).** Needs data payloads and more page weight.
- **C. DuckDB-WASM full browser queries (later).** Needs a hosted data release, a rebuild of stale tables, and hardware that the user wants studied first.

Order: A, then B, then a hardware feasibility study for C.

## Design (phase A)

1. **Taxonomy and colors.** Add `order` to the level list. Each phylum gets a base hue. Subphylum, class and order shades derive from the parent hue. Colors must stay colorblind-safe (reuse the palette checks in `lib/figures.py`). Groups with fewer than `MIN_N` genomes are pooled as "Other".
2. **Feature tables (per genome).** From `asm_stats`, `gene_info`, `gene_proteins`, `gene_introns`, `gene_exons`, `gene_intergenic_distances`, `codon_freq` (61 codons) and `aa_freq` (20 aa). Output small parquet files under a data directory outside git.
3. **Figure modules** (one per topic). Distributions by group at each level, plus gene density, coding fraction, intron and exon counts and lengths, GC, codon usage and aa composition summaries.
4. **Comparisons.**
   - ANOVA or Kruskal-Wallis per statistic by group, with effect size (eta squared or epsilon squared).
   - Regression between statistics (for example genome size against gene count, intron count against gene density). Report slope, R² and n, with phylogeny not controlled. State that caveat on the pages.
   - PCoA on genus representatives (1,671 genera; one genome per genus, picked by best BUSCO then N50; this rule is my call). Distance matrices for codon usage, aa usage and a scaled structure-stat vector. PERMANOVA by phylum and by order. Related views: CA and UMAP.
5. **Figure index.** `figures_index.tsv` with name, topic, level, statistic, n genomes, n groups, effect size.
6. **Selection strategy.** A gallery page lists every figure with filters. A scoring table ranks figures by effect size, coverage and low pooled share. The user ranks from the gallery. Chosen figures move to main pages.
7. **Pages.** Overview, one page per topic, comparisons, ordination, gallery, data and methods. High-level phylum view first, then tabbed breakdown (extend the existing `fig()` helper with an order tab).
8. **Hosting.** A GitHub Actions workflow publishes the rendered `_site` to Pages. A `CNAME` file holds `bfd.fungalgenomes.org`. Figures and small tables are committed, or built by CI from a committed small-table folder. Data larger than a few MB stays out of git. Enabling Pages and the DNS check are the user's steps (repo settings).

## Decisions I made (not yet confirmed)

- Representative rule: best BUSCO, then N50.
- Use PCoA on genus representatives, not on all genomes.
- Mixed order level uses "unassigned" instead of dropping genomes.
- No subagent builder check was run for this document. I checked it myself.

## Left to the builder

- Exact color assignments per phylum.
- Which distance metric for each composition matrix (Bray-Curtis or Euclidean on CLR).
- Page layout details and gallery filter implementation (static HTML or a table widget).

## Open items

- Hardware requirements for phase C: a separate feasibility study. The user's message on this was cut off, so the needs are not yet recorded.
- Confirm that the tables used are current after the 2026-11-16 rebuilds.
