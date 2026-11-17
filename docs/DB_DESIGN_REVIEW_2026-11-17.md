# BFD analytics database: structure and query-performance review (2026-11-17)

**Provenance.** This is the report of a read-only review by a Fable subagent, requested by the project owner on 2026-11-17. The reviewer ran the queries and timings below. The session that asked for the review (Claude Sonnet 5.5) did not re-run them. Treat every number as the reviewer's single-run measurement unless it is marked INFERRED. The report is reproduced as delivered, apart from this header and the status section.

**Why it was requested.** `SELECT is_gpi_anchored, count(*) FROM v_protein_annotation GROUP BY 1` did not finish in 280 s. The owner asked whether better indexes would help and for a review of the database design.

## Status of the recommendations (added 2026-11-17, not part of the review)

| Item | Status on 2026-11-17 |
|---|---|
| Pipeline database `Fungi_BFD_runs/BFD_function/db/BFD.duckdb` has the old `is_gpi_anchored` join | Rebuild submitted as SLURM job 32132750 (`BFD_function/rebuild_pipeline_duckdb.sh`); pending when this was written |
| R1 materialized `protein_annotation` table | Experiment E1 submitted as SLURM job 32132757 (`BFD_function/experiments/E1_protein_annotation/run_E1.sh`); pending when this was written. No change to the pipeline yet |
| Caveats in section 3, items 1 to 7 and 11 | Summarized in `BFD_VIEW_REVIEW.md` section 8 |
| R2 to R10 | Not started |

---

## Test conditions
- Host r11, inside a SLURM job with 4 CPUs and a 16 GB cgroup. It is not an 8-thread node, so every timing below uses 4 threads and an 8 or 9 GB DuckDB memory limit.
- DuckDB CLI 1.1.3 unless marked "1.5.5" (python duckdb 1.5.5).
- Database: `BFD_dashboard/db/BFD.duckdb`, opened read-only. Parquet under `Fungi_BFD_runs/tables/`.
- Timings are single runs on a shared GPFS filesystem. Treat them as approximate.
- Scratch was `$TMPDIR/bfd_review` (node-local). The reviewer deleted it and wrote nothing under `/bigdata`.

Two incidents:
- The reviewer ran `SELECT is_gpi_anchored, count(*) FROM v_protein_annotation WHERE LOCUSTAG='F77FFF1D' GROUP BY 1` as a one-genome test. It is in effect the full join (section 3, item 1). The kernel killed it (exit 137) at the 16 GB cgroup limit, inside a 280 s timeout. It left no temporary files.
- `BFD_dashboard/db/BFD.duckdb.tmp/` is empty and has mtime 19:18 on 2026-11-16. The reviewer's first catalog queries ran at about that time, before it set `temp_directory`. It cannot exclude that they touched the directory mtime.

---

## 1. Would indexes help v_protein_annotation group-bys?

No. The evidence supports the prior stated in the request.

- **No index is possible on the current design.** MEASURED: `duckdb_indexes()` lists only the two indexes on `swissprot_subject`, the one base table. All 59 other objects are views over `read_parquet()`. DuckDB cannot put an ART index on a view or a parquet file.
- **The plan has no step an index could serve.** MEASURED (`EXPLAIN` of the full group-by, 1.1.3 and 1.5.5): 10 `HASH_JOIN` operators, full scans of all 11 inputs, and 5 `HASH_GROUP_BY` operators below the joins.
- **The real cause is that DuckDB does not remove unused LEFT JOINs.** MEASURED from the same plan: the query needs only `gene_proteins` and `predgpi`. The plan still builds hash tables on deeploc (244.8 M rows, three string columns), wolfpsort (244.8 M), the pfam GROUP BY (368 M rows in), idp_summary (146 M), tmhmm, merops, signalp, cazy and targetp. It still reads the unused payload columns (`localizations`, `IDP_fraction`, `PredHel`).
- **The same answer from the two tables it needs takes seconds.** MEASURED:
  ```sql
  SELECT (g.protein_id IS NOT NULL) AS is_gpi_anchored, count(*)
  FROM gene_proteins gp
  LEFT JOIN (SELECT DISTINCT protein_id FROM predgpi WHERE feature='GPI-anchor') g USING(protein_id)
  GROUP BY 1;
  ```
  17.6 s (1.1.3) and 11.3 s (1.5.5). Result: 3,299,703 true, 241,513,215 false.
- **A materialized table makes the group-by sub-second on a 1/16 sample.** MEASURED on LOCUSTAG `F7*` (1,486 genomes, 15,385,240 proteins):

  | Query | Time |
  |---|---|
  | View, `is_gpi_anchored` group-by on the F7 range | 27.6 s |
  | Same group-by on the materialized parquet | 0.09 s |
  | Group by `deeploc_localizations` | 0.24 s |
  | Group by `LOCUSTAG` with 3 aggregates | 0.25 s |
  | 4-column group-by | 0.32 s |
  | `WHERE LOCUSTAG='F77FFF1D'` | 0.012 s |

  The materialized file was on node-local NVMe, not `/bigdata`. INFERRED: linear scaling gives 1.5 to 5 s for all 244.8 M rows. The reviewer did not measure the full table.

Indexes on a copied DuckDB table would help point lookups only. The builder header records this (2026-11-07: one genome in `gene_proteins`, 1.9 s with copies and indexes against 8.5 s with views). They cost a 394 GB database. `in_genome()` now gives 0.05 s for the same lookup with no index (MEASURED).

---

## 2. Ranked recommendations

### R1. Materialize `protein_annotation` as a parquet table at merge time
- **What:** a new MERGE step writes `tables/protein_annotation.parquet`, one row per protein, with LOCUSTAG stored and sorted by protein_id. Keep `v_protein_annotation` as a thin view over it, or as a documented slow path.
- **Benefit:** MEASURED on the 1/16 sample: 27.6 s becomes 0.09 to 0.32 s, about 100 to 300 times faster. Group-by on LOCUSTAG works without the species_prefix trap.
- **Cost:** MEASURED: 132 MB for 15.4 M rows (zstd). INFERRED: about 2.1 GB for all proteins. `protein_id` (48 MB), `length` (30 MB) and `IDP_fraction` (24 MB) dominate; the boolean flags are near 0 MB.
- **Build time:** MEASURED: 27.9 s for the F7 chunk at 9 GB. INFERRED: 16 chunks take about 7.5 minutes at 4 threads.
- **Risk:** the table goes stale if a source table changes without a rebuild. Tie it to the merge manifest. The correctness problems in section 3 get frozen in, so fix them first.
- **Build it in key-range chunks.** MEASURED: the same F7 COPY with `ORDER BY protein_id` failed with "Out of Memory ... 7.4 GiB/7.4 GiB used" at an 8 GB limit. A single unchunked build was not tested. It needs a large node.

### R2. Store the genome key in every table and stop computing it in views
- **What:** write `species_prefix` (or LOCUSTAG) into gene_proteins, gene_transcripts, gene_exons, gene_CDS, gene_introns, gene_trnas and pfam. Remove the `string_split(...)[1]` expression from `mkview`.
- **Benefit:** MEASURED:

  | Query | Time |
  |---|---|
  | `count(*) FROM gene_proteins WHERE species_prefix='F77FFF1D'` (computed) | 11.9 s, 37 s CPU (7.2 s on 1.5.5) |
  | Same with `in_genome(protein_id, ...)` | 0.05 s |
  | `count(*) FROM deeploc WHERE species_prefix='F77FFF1D'` (stored) | 0.5 s |
  | `GROUP BY species_prefix` on gene_proteins (computed) | 13.0 s, 42.8 s CPU |
  | `GROUP BY species_prefix` on wolfpsort (stored, same row count) | 5.9 s, 8.3 s CPU |

- **Cost:** MEASURED: the stored `species_prefix` column is 0.0 MB compressed in every sorted table (dictionary plus RLE).
- **Risk:** `SELECT *` column order changes. The gene_stats builder needs a change and a rerun.

### R3. Make v_protein_annotation safe until R1 exists
- Expose a genome filter that prunes. `WHERE LOCUSTAG = 'X'` on the view filters only the gene_proteins side, after `string_split`. MEASURED (EXPLAIN): the other 10 scans have no filter, so it runs the full-size build sides. `in_genome(protein_id, 'X')` reaches all 11 scans (MEASURED: 2.1 to 4.5 s).
- Options: a table macro `protein_annotation_for(locustag)`, or R2 plus join conditions on both `species_prefix` and `protein_id`. INFERRED that the second option pushes the filter down. Not tested.
- The pfam `GROUP BY protein_id` inside the view is not pushed down or removed. MEASURED: a 368 M row hash aggregate (estimated 184 M groups) runs on every unfiltered query, even when `pfam_domain_count` is not selected. Precompute `pfam_domain_count` per protein at merge time, as part of R1.
- **Cost:** small. **Risk:** low.

### R4. Split wide tables into narrow and wide parts; investigate slow narrow scans on deeploc
- MEASURED: `SELECT species_prefix, count(*) FROM deeploc GROUP BY 1` took 35.8 s, 41.7 s, 38.0 s and 35.2 s (1.5.5). CPU was only 6 to 8 s. The column is 0.4 MB in 1,992 chunks.
- MEASURED: `count(DISTINCT species_prefix)` took 32.4 s on deeploc and 0.24 s on wolfpsort (1.5.5, one session, same row count and row-group count).
- MEASURED: raw `pread` of the 1,992 chunks took 9.6 s on deeploc (4.8 ms per read) and 2.9 s on wolfpsort (1.5 ms per read).
- MEASURED: `GROUP BY localizations` on deeploc took 44.8 s.
- The cause is not determined. INFERRED: per-row-group seek latency on GPFS in a 10 GB file where the 14 `p_*` FLOAT columns hold 9.8 of 10.4 GB.
- **What:** put the label columns (`localizations`, `signals`, `membrane_types`, the two booleans) in R1's table. Keep the `p_*` columns in a separate file. Do the same for `md5checksum`, which is 4.47 of 5.53 GB in gene_proteins and 4.49 of 6.70 GB in gene_transcripts (MEASURED).
- **Benefit:** INFERRED. **Cost:** one more file per split table. **Risk:** queries that need both parts need a join.

### R5. Make sort order explicit in every merge
- MEASURED from parquet footers (row groups whose key range overlaps the previous one):
  - Fully sorted: pfam, deeploc, merops, signalp, tmhmm, swissprot_best, tcdb (0 overlaps); predgpi, wolfpsort (1).
  - Neighbours overlap: gene_proteins 1,036 of 1,977; gene_transcripts 623 of 2,014; gene_CDS 3,104 of 6,083; idp_summary 855 of 1,183; idp 1,252 of 1,897. A one-genome range still opens only 1 or 2 row groups.
  - Not sorted: cazy.overview 75 of 76; cazy.cazymes_hmm 74 of 75; tf_inventory 88 of 89. A one-genome range opens every row group (76 of 76 and 89 of 89).
- `nextflow/sql/merge/` has `ORDER BY` only in pfam.sql, swissprot_*.sql and merge_deeploc.py. wolfpsort, predgpi, tmhmm, merops and signalp are sorted only because of input order.
- **What:** add `ORDER BY species_prefix, protein_id` to every merge SQL. Add a build check on footer overlap.
- **Benefit:** protects the 0.05 to 0.5 s genome lookups. **Cost:** sort time and memory at merge; not measured. **Risk:** low.

### R6. Row-group size: keep the default
- MEASURED: all large tables use about 122,880 rows per group. swissprot uses 1,000,000.
- INFERRED: a genome averages about 10,500 proteins, so one row group holds about 12 genomes. Larger groups would make genome lookups read more and would not speed full scans much. The reviewer has no measurement that supports a change.

### R7. Do not partition by genome
- INFERRED: 23,412 hive partitions per table means 23,412 footers per full scan. The builder comment already records 7 to 10 s per protein for opening the 256 swissprot_hits bucket files through a view. Sorted single files plus range filters already give 0.05 to 0.5 s genome lookups (MEASURED).
- If a partition is wanted, use the first two hex characters (16 values, MEASURED 1,398 to 1,550 genomes each). Not tested.

### R8. Do not adopt a fixed-width integer genome id yet
- MEASURED: 3 of 23,680 LOCUSTAGs have 9 characters (`FF38F669A`, `F3DC5624A`, `FCB937CBA`). Each shares its first 8 characters with another genome (`FF38F669`, `F3DC5624`, `FCB937CB`). A hex-to-UINT32 key would fail or collide.
- **Benefit** is small. MEASURED: the stored VARCHAR key costs 0.0 MB, and `protein_id` costs 190 MB for 245 M rows in sorted files. INFERRED: an integer key would give faster hash joins. Not measured.
- **Cost:** every analysis must map ids back.

### R9. Types and flag columns
- `gene_transcripts.has_start_codon` and `has_stop_codon` are VARCHAR 'TRUE'/'FALSE'/'NULL'. Storage is 7 MB each (MEASURED), so this is a correctness problem, not a speed problem. Write BOOLEAN with real NULL at merge. `v_gene_transcripts` already patches it.
- BIGINT for strand, order and counts, and DOUBLE for scores: INFERRED small benefit from narrower types, because parquet encoding already compacts them (pfam `domain_num` is 50 MB for 368 M rows, MEASURED). Low priority.
- Codec: predgpi, tmhmm, wolfpsort, merops, signalp, cazy, tf_inventory, tcdb and targetP on disk are SNAPPY (MEASURED). The merge SQL now says ZSTD, so these files predate that change. The next merge will recompress them. The size effect was not measured.

### R10. Build workflow
- **Check uniqueness on every table the view joins without DISTINCT.** `check_unique` does not cover deeploc, wolfpsort, tmhmm, targetp, signalp, merops or predgpi. A duplicate protein_id would multiply rows in v_protein_annotation silently. MEASURED on 2026-11-16: 0 duplicate protein_id in wolfpsort (244,790,951 rows), tmhmm (45,707,260), deeploc (244,812,918) and targetp (178,535). The check ran in 16 key ranges, 97 s for all four, under 8 GB.
- **Build one database.** MEASURED: the pipeline copy (`BFD_function/db`, 2026-11-14) has 39 views and still has the unfiltered `FROM predgpi` join, so `is_gpi_anchored` is wrong there. The dashboard copy (2026-11-16) has 59 views and the fix. The pipeline copy's paths go through the `BFD_function/tables` symlink.
- **Set `temp_directory` to `$SCRATCH` for every reader.** By default DuckDB spills to `BFD.duckdb.tmp/` next to the database, which is on `/bigdata`.
- **Use an atomic publish in the pipeline builder.** `build_BFD_duckDB.sh` runs `rm -f "$DB"` first. The dashboard wrapper stages and renames; the pipeline process does not.
- **Add a regression test for the view:** compare `count(*)` of the annotation table with `gene_proteins`, and assert that the fraction of `is_gpi_anchored` is plausible. That test would have caught the GPI bug.
- **Minor:** MEASURED: the database file is 424.5 MiB with 562 of 1,698 blocks free, for one 575 k row table. The reviewer did not find what uses the remaining space.

---

## 3. Wrong or fragile items

1. **`WHERE LOCUSTAG = ...` on v_protein_annotation runs the full join.** MEASURED by EXPLAIN, and by the 16 GB kill above. LOCUSTAG is the column the view offers for this purpose.
2. **The pipeline database copy still has the GPI bug.** See R10.
3. **`is_gpi_anchored = FALSE` also covers "no prediction".** MEASURED (1.5.5, 13.4 s): 17,574,741 proteins have no predgpi row. 2 genomes have none. 1,818 genomes have fewer predgpi rows than proteins. 1,064 genomes lack rows for more than half of their proteins. The largest gap is 100,490 proteins in one genome. In F77FFF1D predgpi has exactly one row per protein (10,110 Chain plus 165 GPI-anchor = 10,275).
4. **`IDP_fraction` is NULL where other flags use 0 or FALSE.** MEASURED: idp_summary is positive-only (F77FFF1D: 6,305 rows for 10,275 proteins, minimum 0.0168, no zero). In the F7 sample 9,195,486 of 15,385,240 proteins have a value. NULL can mean "no disorder" or "not run". The join itself works since the 2026-11-11 rebuild.
5. **`targetp_class = 'OTHER'` is a silent constant.** MEASURED: 11,776 of 15,385,240 proteins in the F7 sample have another value. TargetP covers 125 genomes (from the view review; not recounted).
6. **`tmhmm_helices = 0`, `has_signal_peptide = FALSE` and `wolfpsort_loc` hide missing coverage.** MEASURED: 2 genomes have no wolfpsort rows, and 2,631 genomes differ in wolfpsort row count from gene_proteins (21,967 proteins in total). tmhmm has rows for 23,409 genomes; signalp and predgpi for 23,410; gene_proteins for 23,412.
7. **`left(id, 8)` and 8-character ranges merge two genomes for the three 9-character LOCUSTAGs.** MEASURED: the range `'FF38F669'` to `'FF38F66A'` returns 20,620 proteins; `in_genome(protein_id, 'FF38F669A')` returns 13,799. `in_genome()` and `split_part(id, '_', 1)` are correct. `BFD_VIEW_REVIEW.md` lists `left(id,8)` as a pattern in the analyses.
8. **Sort order is accidental for most tables, and three small tables are unsorted.** See R5.
9. **`mkview` uses `SELECT *` with no schema check** for all tables except the three schema-contract views. A changed parquet schema passes silently.
10. **Views hold absolute paths.** The database works only where `tables/` is at the same path.
11. **Row counts differ slightly from `BFD_VIEW_REVIEW.md`** because of the 2026-11-11 rebuild (gene_proteins 244,812,918 now; 244,988,294 then). Documents that quote the old counts are stale.
12. Items from the earlier review that the reviewer did not recheck: infernal_rfam column shift, the gene_introns off-by-one (since fixed, see `BFD_VIEW_REVIEW.md` section 4.3), the asm_stats column renames.

---

## 4. Experiment plan for the top 3

Use one SLURM job, 16 CPUs, 64 GB, with `SET temp_directory = '$SCRATCH/duckdb_tmp'`. Write test output to `$SCRATCH` first.

### E1 (R1): full materialization, chunked
```sql
SET threads=16; SET memory_limit='48GB'; SET preserve_insertion_order=false;
-- run 16 times, key ranges '' to 'F1', 'F1' to 'F2', ... 'FF' to '~'; one file per chunk
COPY (SELECT * FROM v_protein_annotation WHERE protein_id >= 'F7' AND protein_id < 'F8')
  TO '$SCRATCH/pa/pa_F7.parquet' (FORMAT parquet, COMPRESSION zstd);
```
Then:
```sql
SELECT count(*), count(DISTINCT protein_id) FROM '$SCRATCH/pa/*.parquet';   -- expect 244,812,918 twice
SELECT is_gpi_anchored, count(*) FROM '$SCRATCH/pa/*.parquet' GROUP BY 1;   -- expect 3,299,703 true
SELECT LOCUSTAG, count(*) FROM '$SCRATCH/pa/*.parquet' GROUP BY 1;
SELECT deeploc_localizations, wolfpsort_loc, count(*) FROM '$SCRATCH/pa/*.parquet' GROUP BY ALL;
```
- Record the build time per chunk, the total size, and each query time.
- Expected from the sample: about 28 s per chunk at 4 threads, about 2.1 GB, group-bys in 1.5 to 5 s.
- Then copy one merged, sorted file to `/bigdata` and repeat the queries. The sample timings came from NVMe, so this step measures the GPFS effect.
- Optional: one unchunked `COPY` with a 48 GB limit, to learn whether the full join completes and how long it takes.

### E2 (R2): stored genome key on gene_proteins
```sql
COPY (SELECT split_part(protein_id,'_',1) AS species_prefix, * FROM read_parquet('.../gene_proteins.parquet')
      ORDER BY species_prefix, protein_id)
  TO '$SCRATCH/gene_proteins_sp.parquet' (FORMAT parquet, COMPRESSION zstd);
SELECT count(*) FROM '$SCRATCH/gene_proteins_sp.parquet' WHERE species_prefix='F77FFF1D';   -- now 11.9 s through the view
SELECT species_prefix, count(*) FROM '$SCRATCH/gene_proteins_sp.parquet' GROUP BY 1;        -- now 13.0 s, 42.8 s CPU
```
- Check the footers: 0 overlapping row groups and a species_prefix column near 0 MB.
- The sort of 245 M rows needs the 64 GB node.
- Then rewrite the v_protein_annotation joins to `USING (species_prefix, protein_id)` in a scratch database. Use `EXPLAIN ... WHERE LOCUSTAG='F77FFF1D'` to confirm that the filter reaches all scans.

### E3 (R4): narrow against wide file for deeploc
```sql
COPY (SELECT species_prefix, protein_id, localizations, signals, membrane_types,
             has_signal_peptide, has_mito_transit_peptide
      FROM read_parquet('.../deeploc.parquet') ORDER BY species_prefix, protein_id)
  TO '<test dir on /bigdata>/deeploc_labels.parquet' (FORMAT parquet, COMPRESSION zstd);
SELECT species_prefix, count(*) FROM '<narrow file>' GROUP BY 1;                 -- wide file: 35 to 42 s
SELECT localizations, count(*) FROM '<narrow file>' GROUP BY 1;                  -- wide file: 44.8 s
```
- The narrow file must be on `/bigdata`, because the test is about GPFS behaviour. Run each query twice.
- Also run the wide-file queries with 16 threads. That shows whether the 35 s is per-read latency (it should fall with threads) or something else.
- 16 CPUs and 32 GB are enough.

---

## What the reviewer could not check
- The real runtime and peak memory of the full unfiltered view. It did not rerun it. The only data are the 280 s timeout and its accidental run that was killed at 16 GB.
- Any 8-thread or 16 GB-limit timing. The job had 4 CPUs and a 16 GB cgroup.
- Full-scale materialized group-by times, and any materialized-table timing on `/bigdata`.
- The cause of the slow narrow scans on deeploc.parquet.
- The effect of row-group size, of integer keys on join speed, and of ZSTD against SNAPPY. No test was run.
- Whether `USING (species_prefix, protein_id)` gives pushdown for a LOCUSTAG filter.
- Uniqueness for signalp, merops, predgpi and cazy_overview (the view applies DISTINCT to those), and for the gene_* tables (the build checks those).
- The second half of `BFD_VIEW_REVIEW.md` (lines 271 to 383). The reviewer read lines 1 to 270.
- The merge Python scripts in detail. It checked only their sort and compression clauses.
- MotherDuck. Not queried.

## Files the reviewer read
- `Fungi_BFD/nextflow/bin/build_BFD_duckDB.sh`
- `Fungi_BFD/nextflow/sql/merge/*.sql` (grep for sort and compression only)
- `Fungi_BFD/nextflow/bin/merge_deeploc.py` (grep only)
- `BFD_dashboard/workflows/build_bfd_db.sh`
- `BFD_dashboard/sql/bfd_extra_views.sql` and `bfd_compat_views.sql`
- `BFD_dashboard/docs/BFD_VIEW_REVIEW.md` (lines 1 to 270)
- `BFD_dashboard/db/BFD.duckdb` and `Fungi_BFD_runs/BFD_function/db/BFD.duckdb` (read-only)
- Parquet footers in `Fungi_BFD_runs/tables/`
