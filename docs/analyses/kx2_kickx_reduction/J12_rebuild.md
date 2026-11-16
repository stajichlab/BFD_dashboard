# J12: end-to-end rebuild of Analysis E (H1-H4) from the local BFD database

Date: 2026-11-15. Status: **rebuild ran; headline H1 and H2 numbers reproduce, most other tables differ slightly, and two input differences were found.**
The earlier Analysis E numbers had never been reproduced from empty scratch (only later steps were re-run). J12 does that on the HPCC, from the saved SQL.

## What was run

| Item | Value |
|---|---|
| Job | `workflows/atlas/J12_kx2_rebuild.sh`, SLURM job 32074412, 16 cores, 3 min 45 s |
| Driver | `bin/atlas/j12_kx2_rebuild.py` (runs the saved SQL in `sql_upstream/` and `sql_h1..h4/` in order; fixture test `bin/atlas/test_j12_kx2_rebuild.py`: 4 passed on the HPCC, job 32102350) |
| Input | the shared local DuckDB `BFD.duckdb` (read-only), the local replacement for MotherDuck introduced in `workflows/build_bfd_db.sh` |
| Scratch | empty DuckDB file in node-local `$SCRATCH`, removed when the job ends; only small tables are written to the project folder |
| Outputs here | [`j12_rebuild/`](j12_rebuild/): `run_log.tsv` (every statement and its time), `table_counts.tsv`, `table_fingerprint.tsv` (per column: non-null count and sum), `results/*.tsv` (result queries of the upstream, H1 and H2 stages) |

The saved SQL did not include the definition of the genome set `kx_genomes` or the species/assembly staging tables (`j12_species`, `j12_asm_stats`). `sql_upstream/00-02` **reconstruct** them; the reconstruction is the main source of difference below.

## Row counts: 27 of 60 tables identical, 33 differ

Identical in row count to the MotherDuck originals: `kw_pfam_matrix` (32,265,977), `kw_rep_genomes`, `kx2_h1_decomp`, `kx2_h1_fam_gen`, `kx2_h1_genus`, `kx2_h1_arch_metrics`, `kx2_h1_fusion`, `kx2_h1_genomes`, `kx2_h2_delta`, `kx2_h2_lin`, `kx2_h3_codmap`, `kx2_h3_cliff`, `kx2_h3_grpmed`, most of the H4 tables (`kx2_h4_loss`, `kx2_h4_overlap`, `kx2_h4_gcls`, `kx2_h4_resid`, `kx2_h4_gp_all`, `kx2_h4_lostset`, `kx2_h4_pfcat`, `kx2_h4_tf`, ...) and others (full list: compare `j12_rebuild/table_counts.tsv`).

| table | MotherDuck | J12 rebuild | change |
|---|---|---|---|
| `kx_genomes`, `kx2_h2_genomes`, `kx2_h3_aa`, `kx2_h3_cod`, `kx2_h3_genome`, `kx2_h3_trna` | 1,919 / 1,919 / 1,926 each | 1,945 / 1,945 / 1,952 each | +1.35% |
| `kx2_h1_genome_all`, `kx2_h1_genome_metrics`, `kx2_h1_longfam` | 740 | 728 | -1.62% |
| `kx2_h1_prot` | 4,998,294 | 4,954,966 | -0.87% |
| `kx2_h1_cliff` | 300 | 280 | -6.67% |
| `kx2_h2_tx`, `kx2_h2_ig` | 12,423,748; 11,686,057 | 12,479,951; 11,760,258 | +0.45%; +0.63% |
| `kx2_h4_meta`, `kx2_h4_nprot` | 7,411 | 7,393 | -0.24% |
| `kx2_h4_rep_gen` | 4,064 | 4,052 | -0.30% |
| other H2/H3/H4 tables | within 0.0-1.4% | | see `table_counts.tsv` |

## The headline numbers (observed)

| Claim in `FUNCTIONAL_ANALYSES.md` section E | Rebuild | Verdict |
|---|---|---|
| H2: core Kickxellales vs Mucorales, deficit 15 Mb, about 70% gene number, about 22% intergenic, 9% introns | 15.0 Mb; 70%, 22%, 9% (`results/h2__06_decomposition__q1.tsv`, quality-filtered row `Q`) | **reproduced exactly** |
| H2: Saccharomycotina vs Pezizomycotina 27.5 Mb, introns 2% | 27.5 Mb; 72% gene number, 23% intergenic, 2% introns | **reproduced** |
| H1: mean protein length unchanged in core Kickxellales (delta +0.09) | `mean_len` +0.09 | **reproduced** |
| H1: multidomain fraction, core Kickxellales vs Mucorales delta -0.58 | `f_multi` -0.58; Harpellales -0.64, Dimargaritales -0.33, Microsporidia -1.00 | **reproduced** |
| H1: -0.30 "on non-partial genes" | `f_multi_np` is identical to `f_multi` (-0.58) in every group | **not reproduced; cause below (partial flag)** |
| H4: copy loss in Kickxellales/yeasts, overlap enrichment (41% / 59%, 0.62-0.89 observed over expected, overlap 1.87 with Harpellales) | table-level sums only (see next paragraph); the J12 result queries for H3/H4 were not written out as separate files | **consistent, not checked number by number** |
| H3: ENC residual +0.054, tRNA match 0.47 vs 0.61, GC3 67% | table-level sums only | **consistent, not checked number by number** |

Table-level checks for H3/H4 (column sums in `table_fingerprint.tsv` against the same sums computed on the MotherDuck tables): `kx2_h4_loss` has the same `n_fam` (188,872) and `base_copies` (557,308); `loss_copies` 250,561 vs 250,653 (-0.04%) and `complete_copies` 88,337.5 vs 88,366 (-0.03%).
`kx2_h4_overlap` has the same enrichment sum (11.62) and `lostB` (7,481), with observed overlap 2,282 vs 2,306 (-1.0%). `kx2_h3_grpmed` has the same GC3 sum (1,329.2), with `enc_dep` 2.014 vs 2.005 and `r_supply` 10.67 vs 10.45. `kx2_h1_cliff_ref` cliff sum -0.73 vs -0.66.
These are sums over a few dozen to a few hundred rows, so they show the tables are close, not that every cell matches.

## Why the numbers differ (observed, and what is still unknown)

1. **The genome set was reconstructed, not saved.** The original `kx_genomes` definition was only described in prose. The rebuilt table matches the MotherDuck counts in every group relevant to the Kickxellomycota questions (Kickxellomycota groups, Mucorales, Zoopagales, Saccharomycotina classes, Pezizomycotina classes except one), but has **26 more genomes** in four unrelated groups:
   Microsporidia 102 vs 83 (+19), chytrids + Monoblepharidomycota 63 vs 60 (+3), *Malassezia* 24 vs 21 (+3), Pezizomycotina Sordariomycetes 280 vs 279 (+1). The reason (annotation-status or species-deduplication rule) was not found. Microsporidia enter Analysis E only as the one lineage with whole-family simplification, so this affects H1 and H4 for them, not the Kickxellomycota comparisons.
2. **`gene_transcripts.is_partial` is very different in the two databases.** In MotherDuck `bfd`, 36,237,955 of 249,546,610 transcripts (14.5%) are partial, and the flag is exactly "missing start or stop codon". In the local `BFD.duckdb`, **357,764 of 249,366,496 (0.14%)** are partial (again identical to missing start or stop). The rebuild's non-partial metrics therefore equal the all-protein metrics, which explains the non-reproduced `f_multi_np` and the small `mean_len_np` shifts, and it also changes the H1 partial-gene
   fractions. The local copy also has 180,114 fewer transcripts (249,366,496 vs 249,546,610), so the two databases are different snapshots or builds. Which copy has the right partial flags was not determined; `gene_transcripts` is loaded from the shared parquet, not built by `workflows/build_bfd_db.sh`. This should go to the owner of the BFD build as an issue.
3. **Not separated:** how much of each remaining shift comes from the extra genomes (1) versus the transcript snapshot (2).

## Interpretation (not shown by the data)

The structure of Analysis E reproduces from saved SQL: the H2 accounting and the H1 architecture results are stable to a different database build and a reconstructed genome set. The one claim that depends on the partial-gene flag (-0.30 on non-partial genes) cannot be checked until the flag question is settled; section E should keep it, marked as unreproduced.

## Limits

- Only H1 and H2 results have per-query output files; H3 and H4 were compared through table-level sums (and one-sided: sums close does not mean each row equal).
- The reconstructed `kx_genomes` was validated against group counts and a per-group LOCUSTAG hash, not by an independent definition.
- No statistical tests: the comparison is a replication of deterministic SQL, so any difference is a difference in input data or rules.

## Next

1. Open an issue for the BFD build owner: `gene_transcripts.is_partial` (and `has_start_codon`/`has_stop_codon`) differ between MotherDuck and the local `BFD.duckdb` (14.5% vs 0.14% partial).
2. Find why 26 more genomes enter the rebuilt `kx_genomes`, then freeze the definition in the repo.
3. Write the H3 and H4 result queries to files in J12, so the headline numbers can be compared directly, and re-run after the partial-flag question is answered.
