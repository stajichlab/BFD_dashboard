# bfd_scratch: scratch tables on the local database

`bfd_scratch` holds intermediate tables from exploratory analyses. It replaces the MotherDuck
database of the same name. The BFD data itself is never changed.

## Layout

| Catalog | File | Access | Content |
|---|---|---|---|
| `bfd` | `db/BFD.duckdb` | read-only | VIEWs over `Fungi_BFD_runs/tables/*.parquet` and the SwissProt macros |
| `bfd_scratch` | `db/bfd_scratch.duckdb` | read-write | analysis tables, one prefix per analysis |

Both files are ignored by git. `workflows/build_bfd_db.sh` makes `BFD.duckdb`.
`bin/bfd_scratch.py init` makes the scratch file.

The catalog names match MotherDuck. SQL such as `CREATE TABLE bfd_scratch.kw_rep_genomes AS ...`
runs unchanged. Bare names such as `pfam` resolve to `bfd`.

## Use

```python
import sys; sys.path.insert(0, "/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/lib")
import bfd_connect
con = bfd_connect.connect(threads=16, memory="60GB", tmp="/bigdata/.../tmp")  # inside a SLURM job
con.execute("CREATE TABLE bfd_scratch.kw_rep_genomes AS SELECT LOCUSTAG FROM species")
```

```bash
/usr/bin/python3.12 bin/bfd_scratch.py init
/usr/bin/python3.12 bin/bfd_scratch.py list
/usr/bin/python3.12 bin/bfd_scratch.py sql "SELECT count(*) FROM bfd_scratch.kw_rep_genomes"
/usr/bin/python3.12 bin/bfd_scratch.py drop kw_
```

## Conventions

1. Name each table `<prefix>_<name>`. The prefix is the analysis (`kw_` kingdom domains,
   `kx_` Kickxellomycota). `bfd_connect.check_prefix()` tests a name.
2. Write nothing into `bfd`. It is attached read-only.
3. DuckDB allows one writer per file. Two jobs that write to the same scratch file will fail
   with a lock error. For parallel jobs, set `BFD_SCRATCH_DB` to one file per job. Many
   readers are safe with `scratch_read_only=True`.
4. A scratch table is a cache, not a result. Keep the SQL that builds it in the analysis
   folder (`docs/analyses/<name>/`). Export results that must last to parquet or TSV.
5. Drop tables when the analysis ends (`bfd_scratch.py drop <prefix>_`).

## Differences from MotherDuck to check

- View names come from the pipeline builder and can differ from the MotherDuck table names.
  The pipeline uses `cazy`, `signalp`, `targetp`, `aa_frequency`. The MotherDuck loader used
  the parquet stems `cazy_overview`, `signalp_signal_peptide`, `targetP`, `aa_freq`. Compare
  `duckdb_views()` in the new `BFD.duckdb` with the table names used in
  `docs/analyses/*/README.md` after the first build. This has not been checked yet.
- MotherDuck `bfd` held `COMMENT ON TABLE` descriptions (`bin/md_load_bfd.py`,
  `DEFAULT_COMMENTS`). The pipeline views have none.
- Existing `bfd_scratch` tables on MotherDuck are not copied. Rebuild them from the SQL in
  the READMEs, or export them with `COPY ... TO 'x.parquet'` while the account is open.
