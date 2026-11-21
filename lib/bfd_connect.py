"""Open the local BFD database plus the writable bfd_scratch database.

This replaces the MotherDuck setup (database ``bfd`` + database ``bfd_scratch``).
SQL written for MotherDuck keeps working because the catalog names are the same:

    bfd.pfam                      read-only view over tables/*.parquet
    bfd_scratch.kw_rep_genomes    writable scratch table

Rules (see docs/BFD_SCRATCH.md):
  * ``bfd`` is attached READ_ONLY. Nothing here can change it.
  * Scratch tables must start with an analysis prefix (``kw_``, ``kx_`` ...).
  * DuckDB allows one writer per database file. For parallel SLURM jobs, give each
    job its own scratch file with BFD_SCRATCH_DB or the ``scratch_db`` argument.

Paths come from arguments, else the environment, else the defaults below:
    BFD_DB          path to BFD.duckdb         (output of workflows/build_bfd_db.sh)
    BFD_SCRATCH_DB  path to the scratch file   (created on first use)
    BFD_LABELS_DB   path to BFD_labels.duckdb  (optional; built by workflows/build_labels_db.sh)

The labels file holds the dashboard's analysis layer (v_genome_labels, v_genome_traits, ref_* tables). When given, it is
attached READ_ONLY as ``labels`` and searched before ``bfd``, so unqualified ``v_genome_labels`` resolves to it.
Without it nothing changes. The file records which pipeline database it was built from; connect() warns when
that database is newer.

Nothing runs at import time.
"""
from __future__ import annotations

import os
import re
import warnings
from pathlib import Path

import duckdb

_ROOT = Path("/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard")
DEFAULT_BFD_DB = _ROOT / "db" / "BFD.duckdb"
DEFAULT_SCRATCH_DB = _ROOT / "db" / "bfd_scratch.duckdb"

_PREFIX_RE = re.compile(r"^[a-z][a-z0-9]{0,11}_")


def connect(
    bfd_db: Path | str | None = None,
    scratch_db: Path | str | None = None,
    *,
    scratch: bool = True,
    scratch_read_only: bool = False,
    labels_db: Path | str | None = None,
    threads: int | None = None,
    memory: str | None = None,
    tmp: Path | str | None = None,
) -> duckdb.DuckDBPyConnection:
    """Return an in-memory connection with ``bfd`` (read-only) and ``bfd_scratch`` attached.

    scratch=False           attach ``bfd`` only.
    scratch_read_only=True  attach scratch read-only (safe for many concurrent readers).
    labels_db               attach this labels file read-only as ``labels`` (else BFD_LABELS_DB, else none).
    threads / memory / tmp  set DuckDB limits. Use tmp on /bigdata or $SCRATCH inside a job.
    """
    bfd_path = Path(bfd_db or os.environ.get("BFD_DB") or DEFAULT_BFD_DB)
    if not bfd_path.exists():
        raise FileNotFoundError(
            f"{bfd_path} not found. Build it with: sbatch workflows/build_bfd_db.sh "
            "(or set BFD_DB to the pipeline copy under Fungi_BFD_runs/BFD_function/db/)."
        )
    con = duckdb.connect()  # in memory; the real data is in the attached files
    if threads:
        con.execute(f"SET threads={int(threads)}")
    if memory:
        con.execute(f"SET memory_limit='{memory}'")
    if tmp:
        con.execute(f"SET temp_directory='{tmp}'")
    con.execute(f"ATTACH '{bfd_path}' AS bfd (READ_ONLY)")
    if scratch:
        sp = Path(scratch_db or os.environ.get("BFD_SCRATCH_DB") or DEFAULT_SCRATCH_DB)
        if scratch_read_only:
            if not sp.exists():
                raise FileNotFoundError(f"{sp} not found (scratch_read_only needs an existing file)")
            con.execute(f"ATTACH '{sp}' AS bfd_scratch (READ_ONLY)")
        else:
            sp.parent.mkdir(parents=True, exist_ok=True)
            con.execute(f"ATTACH '{sp}' AS bfd_scratch")
    con.execute("USE bfd")  # unqualified names (pfam, species) resolve to bfd, as before
    labels = labels_db or os.environ.get("BFD_LABELS_DB")
    if labels:
        lp = Path(labels)
        if not lp.exists():
            raise FileNotFoundError(f"{lp} not found. Build it with: sbatch workflows/build_labels_db.sh")
        con.execute(f"ATTACH '{lp}' AS labels (READ_ONLY)")
        con.execute("SET search_path = 'labels.main,bfd.main'")
        _warn_if_stale(con, bfd_path, lp)
    return con


def _warn_if_stale(con: duckdb.DuckDBPyConnection, bfd_path: Path, labels_path: Path) -> None:
    """Warn when the labels were built from a pipeline file older than the one attached now."""
    try:
        built_from = con.execute("SELECT source_mtime_epoch FROM labels.labels_meta").fetchone()[0]
    except duckdb.Error:
        warnings.warn(f"{labels_path} has no labels_meta table; cannot tell which pipeline database it came from", stacklevel=3)
        return
    if built_from is not None and int(bfd_path.stat().st_mtime) > int(built_from):
        warnings.warn(
            f"{labels_path} was built from a pipeline database older than {bfd_path}; "
            "rebuild it with: sbatch workflows/build_labels_db.sh", stacklevel=3)


def check_prefix(table: str) -> str:
    """Return the analysis prefix of a scratch table name, or raise ValueError."""
    m = _PREFIX_RE.match(table)
    if not m:
        raise ValueError(
            f"scratch table '{table}' needs an analysis prefix, e.g. 'kw_{table}' "
            "(lowercase letters/digits, up to 12 chars, then underscore)"
        )
    return m.group(0)


def scratch_tables(con: duckdb.DuckDBPyConnection) -> list[tuple[str, str, int]]:
    """List (table, prefix, estimated_rows) in bfd_scratch."""
    rows = con.execute(
        "SELECT table_name, estimated_size FROM duckdb_tables() "
        "WHERE database_name = 'bfd_scratch' ORDER BY table_name"
    ).fetchall()
    out = []
    for name, est in rows:
        m = _PREFIX_RE.match(name)
        out.append((name, m.group(0) if m else "(no prefix)", est))
    return out
