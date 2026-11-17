#!/usr/bin/bash -l
#SBATCH -p epyc -c 2 --mem 24G --time 3:00:00 -J qs_figures
#SBATCH -o %x.%j.log
# Build all summary figures and tables for the Quarto site from the per-genome parquet tables.
# Submit from the repo root: sbatch workflows/summary_figures.sh   (override ROOT, DATA, OUT, PY with env vars)
# The SLURM log (qs_figures.<jobid>.log) is written in the submit directory.
set -euo pipefail
ROOT=${ROOT:-${SLURM_SUBMIT_DIR:?submit from the repo root or set ROOT}}
DATA=${DATA:-/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/summary_tables}
OUT=${OUT:-/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/summary_figures_real}
PY=${PY:-/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/atlas_env/venv/bin/python}
: "${SCRATCH:?SCRATCH not set}"
export MPLCONFIGDIR="$SCRATCH/mplconfig"
mkdir -p "$OUT" "$MPLCONFIGDIR"
cd "$ROOT"
START=$(date +%s)
$PY bin/summary/make_figures.py "$DATA" "$OUT"
echo "elapsed_seconds=$(( $(date +%s) - START ))"
