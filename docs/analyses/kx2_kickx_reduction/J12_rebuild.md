# J12: end-to-end rebuild of Analysis E (H1-H4) from the local BFD database

Date: 2026-11-15. Status: **rebuild ran; headline H1 and H2 numbers reproduce, most other tables differ slightly; one input difference (a reconstructed genome set) and one fixed bug in the old copy (the partial-transcript flag) explain the rest.** MotherDuck is being phased out, so the local `BFD.duckdb` is the reference and the MotherDuck tables below are the historical originals.
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
| H1: -0.30 "on non-partial genes" | `f_multi_np` is identical to `f_multi` (-0.58) in every group | **superseded**: the original used a partial flag with a case bug (below); with the corrected flag, non-partial and all proteins give the same result |
| H4: copy loss in Kickxellales/yeasts, overlap enrichment (41% / 59%, 0.62-0.89 observed over expected, overlap 1.87 with Harpellales) | table-level sums only (see next paragraph); the J12 result queries for H3/H4 were not written out as separate files | **consistent, not checked number by number** |
| H3: ENC residual +0.054, tRNA match 0.47 vs 0.61, GC3 67% | table-level sums only | **consistent, not checked number by number** |

Table-level checks for H3/H4 (column sums in `table_fingerprint.tsv` against the same sums computed on the MotherDuck tables): `kx2_h4_loss` has the same `n_fam` (188,872) and `base_copies` (557,308); `loss_copies` 250,561 vs 250,653 (-0.04%) and `complete_copies` 88,337.5 vs 88,366 (-0.03%).
`kx2_h4_overlap` has the same enrichment sum (11.62) and `lostB` (7,481), with observed overlap 2,282 vs 2,306 (-1.0%). `kx2_h3_grpmed` has the same GC3 sum (1,329.2), with `enc_dep` 2.014 vs 2.005 and `r_supply` 10.67 vs 10.45. `kx2_h1_cliff_ref` cliff sum -0.73 vs -0.66.
These are sums over a few dozen to a few hundred rows, so they show the tables are close, not that every cell matches.

## Why the numbers differ (observed, and what is still unknown)

1. **The genome set was reconstructed, not saved.** The original `kx_genomes` definition was only described in prose. The rebuilt table matches the MotherDuck counts in every group relevant to the Kickxellomycota questions (Kickxellomycota groups, Mucorales, Zoopagales, Saccharomycotina classes, Pezizomycotina classes except one), but has **26 more genomes** in four unrelated groups:
   Microsporidia 102 vs 83 (+19), chytrids + Monoblepharidomycota 63 vs 60 (+3), *Malassezia* 24 vs 21 (+3), Pezizomycotina Sordariomycetes 280 vs 279 (+1). The reason (annotation-status or species-deduplication rule) was not found. Microsporidia enter Analysis E only as the one lineage with whole-family simplification, so this affects H1 and H4 for them, not the Kickxellomycota comparisons.
2. **The partial-transcript flag differs because MotherDuck had a case bug, which was fixed before the local database was built (resolved 2026-11-15).** In MotherDuck `bfd`, 36,237,955 of 249,546,610 transcripts (14.5%) are partial, in almost every genome (16,262 of 23,452 genomes at 10-50% partial, one genome at or below 2%), and the flag is exactly "missing start or stop codon". In the shared parquet and the local `BFD.duckdb`, 357,764 of 249,366,496 (0.14%) are partial
   (median per-genome partial fraction 0.0001; 23,322 of 23,412 genomes below 5%, 30 at 10% or more, maximum 38%). The cause: `nextflow/bin/build_genestats_table.py` compared the first CDS codon with `"ATG"` and the last with the stop codons on text read from the genome, so any soft-masked (lower-case) start or stop made the transcript "partial". Commit `fb1d7cd` in `stajichlab/Fungi_BFD` (2026-11-11 12:00, now on `main`) reads the genome upper-case
   (`Fasta(dna, sequence_always_upper=True)`). All 23,412 per-genome `gene_transcripts` and `gene_introns` caches were regenerated on 2026-11-11 between 15:02 and 18:31 (none older than the fix), the merged `gene_*` parquet files are dated 23:26 and `BFD.duckdb` 23:41; all `gene_introns` splice sites are upper-case (508,436,782 rows, none lower-case). **The local database is the corrected one and the MotherDuck copy is stale.** The rebuild's non-partial metrics therefore equal the all-protein metrics
   (the earlier "non-partial" restriction had removed the 14.5% of proteins that were only "partial" because of lower-case start or stop codons). The local copy also has 180,114 fewer transcripts than MotherDuck (249,366,496 vs 249,546,610) and 23,412 genomes vs 23,452, so the two were also different snapshots.
3. **Not separated:** how much of each remaining shift comes from the extra genomes (1) versus the transcript snapshot (2).

## Interpretation (not shown by the data)

The structure of Analysis E reproduces from saved SQL: the H2 accounting and the H1 architecture results are stable to a different database build and a reconstructed genome set. The one claim that depended on the partial-gene flag, "-0.30 on non-partial genes", and the H1 correlation of the multidomain fraction with the partial-gene fraction (r = -0.70) were computed with the buggy flag, so they are **superseded**, not just unreproduced: they should be recomputed on the corrected flag (the first is already answered by the rebuild: -0.58, the same as all proteins).

## Limits

- Only H1 and H2 results have per-query output files; H3 and H4 were compared through table-level sums (and one-sided: sums close does not mean each row equal).
- The reconstructed `kx_genomes` was validated against group counts and a per-group LOCUSTAG hash, not by an independent definition.
- No statistical tests: the comparison is a replication of deterministic SQL, so any difference is a difference in input data or rules.

## Next

1. Recompute the H1 correlation with the partial-gene fraction (r = -0.70 in section E) on the corrected flag; the 14.5% partial rate it was computed on was a case artifact. Check the 30 genomes that still have 10% or more partial transcripts (maximum 38%); they are probably real (fragmentary gene models), but this was not tested.
2. Find why 26 more genomes enter the rebuilt `kx_genomes`, then freeze the definition in the repo.
3. Write the H3 and H4 result queries to files in J12, so the headline numbers can be compared directly, and re-run after the partial-flag question is answered.
