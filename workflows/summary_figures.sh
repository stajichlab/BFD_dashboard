#!/usr/bin/bash -l
#SBATCH -p epyc -c 8 --mem 48G --time 3:00:00 -J qs_figures
#SBATCH -o /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/summary_figures_real/qs_figures.%j.log
# Build all summary figures and tables for the Quarto site from the per-genome parquet tables.
# Submit: sbatch workflows/summary_figures.sh   (override ROOT, DATA, OUT, PY with env vars)
set -euo pipefail
ROOT=${ROOT:-/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/.claude/worktrees/quarto-site}
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
