#!/usr/bin/bash -l
#SBATCH -p epyc -c 16 --mem 64G --time 4:00:00 -J qs_tables
#SBATCH -o %x.%j.log
# Build the per-genome tables for the Quarto summary from the local BFD.duckdb.
# Submit from the repo root: sbatch workflows/summary_tables.sh
# The SLURM log (qs_tables.<jobid>.log) is written in the submit directory.
# Heavy gene-table group-bys run here, not on the login node. DuckDB temp goes to node-local $SCRATCH.
set -euo pipefail
ROOT=${ROOT:-${SLURM_SUBMIT_DIR:?submit from the repo root or set ROOT}}
OUT=${OUT:-/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/summary_tables}
PY=${PY:-/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/atlas_env/venv/bin/python}
: "${SCRATCH:?SCRATCH not set}"
mkdir -p "$OUT"
cd "$ROOT"
$PY bin/summary/build_summary_tables.py "$OUT" --threads "${SLURM_CPUS_PER_TASK:-16}" --memory 56GB
ls -l "$OUT"
