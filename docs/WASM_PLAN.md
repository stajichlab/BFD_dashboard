# BFD Dashboard: Review and Plan for a WASM + LLM Plotting Interface

Date: 2026-11-08

## 1. Summary

The repo has two generations of tools:

- **Generation 1 (Python, static):** `bin/make_dashboard.py` reads `BFD.duckdb`, computes ordinations and statistics for a taxon subset, and writes one offline `explorer.html`.
- **Generation 2 (browser, started):** `bin/rewrite_parquet.py` prepares parquet files for HTTP range requests. `tests/wasm_smoke.html` benchmarks DuckDB-WASM queries. `web/chat.html` lets an LLM write SQL that DuckDB-WASM runs in the browser.

The plan finishes Generation 2. It has four layers:

1. A versioned data release built on HPCC.
2. Analysis of subsets in the browser.
3. Plots that the LLM writes as Vega-Lite specs.
4. One app that holds the explorer views and the chat.

## 2. What each tool does now

| Tool | Purpose | Status |
|---|---|---|
| `bin/make_dashboard.py` + `lib/` | Selects a taxon subset. Builds feature matrices (assembly, codon, amino acid, Pfam, CAZy, MEROPS, localization). Runs PCoA or CA, optional NMDS, t-SNE and UMAP, and PERMANOVA. Builds differential domain tables. Writes `cache/payload.json` and `explorer.html`. | Works for small subsets. The code assumes about 11–40 genomes. |
| `templates/explorer_template.py` | Draws the overview, ordination, composition, differential and phylogeny views in plain JS and SVG. Uses no CDN. | Works. All values are computed in advance. |
| `workflows/*.sh` | Example SLURM scripts for the upstream pipelines (BFD, NovInvenio, phyling). | Not dashboard code. `run_BFD.sh` is specific to one dataset (11 Exophiala genomes). |
| `bin/rewrite_parquet.py` | Re-sorts large tables by genome or by domain with ZSTD. Builds small per-genome summary tables. | Covers only genome, Pfam, protein length and intergenic data. |
| `tests/wasm_smoke.html` | Runs 8 benchmark queries in DuckDB-WASM. | No timings are saved in the repo. |
| `web/chat.html` | LLM writes SQL. DuckDB-WASM runs it. Chart.js draws bar, line or scatter charts. Allows only SELECT, checks with EXPLAIN, makes one repair attempt. | Working prototype. |

## 3. Problems found

1. **The Python builder cannot run on the full database.** The `species` table has 23,680 genomes. A dense distance matrix at that size needs 23,680² × 8 bytes ≈ 4.5 GB. Eigendecomposition is O(n³). This is a calculation, not a test run. A full-database build from `config/BFD_dashboard.yaml` is not practical.
2. **The demo build is thin.** The committed `cache/payload.json` is a 30-genome Fusarium subset with only the `assembly` feature set.
3. **The WASM data covers fewer tables than the Python dashboard.** `rewrite_parquet.py` leaves out:
   - codon and amino-acid frequency
   - CAZy and MEROPS
   - localization (SignalP, TargetP, WoLF PSORT, PredGPI, IDP)
   - BUSCO lineage, telomeres, WGD Ks, TF inventory, SwissProt

   `aa_freq.parquet` covers 23,452 genomes and is only 2 MB, but the browser cannot query it.
4. **The coverage numbers in the chat prompt may be out of date.** `chat.html` has them typed in (Pfam about 745 genomes, protein about 770, intergenic about 7,800). The snapshot tables were rebuilt on Nov 5–7, after `chat.html` last changed. Pfam coverage was not re-measured.
5. **The model never sees the query result.** It writes the explanation and chart spec before the query runs. It cannot notice an empty or wrong result.
6. **The chart spec is too small.** It has only `{type, x, y}`. There is no color by group, log scale, histogram, box plot, heatmap or facet.
7. **The query timeout may not work.** `conn.cancelSent()` is called for a query started with `conn.query()`. This was not checked. Test it.
8. **Housekeeping:** the untracked file `,` is an exact copy of `workflows/run_BFD.sh`.

## 4. Target design

### Layer 1: data release (Python, on HPCC)

- Extend `rewrite_parquet.py` to write a versioned release directory.
- Add small per-genome tables:
  - codon and amino-acid frequency
  - CAZy and MEROPS family counts
  - localization fractions
  - BUSCO
- Write `manifest.json` with, for each table: row count, genome coverage, schema. Also record build date and source database.
- Compute on HPCC what a browser cannot:
  - whole-database embeddings (for example PCoA or UMAP of Pfam profiles over all genomes), stored as parquet
  - the phylogeny as Newick plus a tip table
- `make_dashboard.py` then becomes only this release step.

### Layer 2: analysis in the browser, per subset

- `lib/ordination.py` needs only numpy, pandas and scipy. Pyodide includes all three.
- The same tested PCoA, CA and PERMANOVA code can run in a Web Worker on a subset that the user selects with SQL.
- Alternative: a JS port of about 150 lines. It loads faster, but the code then exists twice.
- Recommendation: try Pyodide first. Measure:
  - Pyodide load time
  - run time at 100, 500 and 1,000 genomes
- No measurements exist yet for either option.

### Layer 3: plots written by the LLM

- Replace the `{type, x, y}` spec with **Vega-Lite**.
  - It is declarative JSON. The page runs no model-written code.
  - It can be checked against a schema.
  - It supports color, facets, log scales, histograms, box plots and heatmaps.
- Remove any `data.url` from the spec. Feed the DuckDB result in as inline data.
- Change the single JSON reply into a tool-use loop. Tools:
  - `run_sql`: returns columns, row count, summary statistics and the first N rows. The full rows stay in the browser.
  - `describe_table`
  - `plot(vega_lite)`
  - `ordinate(subset_sql, feature_set, method)`: calls Layer 2.
- The model then sees results before it explains them.
- Each chart shows its sample size *n* and table coverage, taken from the manifest.

### Layer 4: one app

- Rebuild the current explorer tabs as fixed queries plus fixed Vega-Lite specs in the same WASM page.
- The subset picker becomes a SQL `WHERE` clause on the `genome` table.
- The chat sits next to the views.

## 5. Phases

| Phase | Work | Check before moving on |
|---|---|---|
| 1. Clean up | Delete `,`. Add a demo dataset (from `TODO.md`). Choose hosting. | Demo builds end to end. |
| 2. Data release v1 | Manifest and missing summary tables. | Run `wasm_smoke.html` against the release and save the timings. |
| 3. LLM plotting | Vega-Lite and the tool loop in `chat.html`. | A fixed set of test questions gives valid specs and correct numbers. |
| 4. Subset analysis | Ordination and PERMANOVA in Pyodide. | Browser results match the Python results on the same subset. Size limits are measured. |
| 5. One app | Port the explorer views. Add whole-database embeddings. | Views match the old `explorer.html` on the demo subset. |

## 6. Open questions (no data yet)

- **Hosting.** The full tables are over 100 GB (Pfam alone is 20 GB). They need a static host with HTTP range requests and CORS. The code says that browser (CORS) access to the nrp.ai endpoint is not confirmed.
- **Query speed.** Speed of DuckDB-WASM range reads on the 20 GB Pfam tables is unknown. Memory use of large joins is unknown. The wasm32 memory limit is 4 GB.
- **API key.** Option A: each user brings their own key, as the page does now. Option B: run a proxy.
