#!/usr/bin/bash -l
#SBATCH -p epyc -c 4 --mem 32G --time 3:00:00 -J summary_ecology
#SBATCH -o %x.%j.log
# Ecology (FungalTraits lifestyle and growth form) in the summary ordinations: global and per-clade PCoA and PERMANOVA.
# Run after workflows/summary_figures.sh. Needs a BFD.duckdb built with trait columns (workflows/build_bfd_db.sh).
# Submit from the repo root: sbatch workflows/summary_ecology.sh   (override ROOT, DATA, DB, OUT, PY with env vars)
set -euo pipefail
ROOT=${ROOT:-${SLURM_SUBMIT_DIR:?submit from the repo root or set ROOT}}
DATA=${DATA:-/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/summary_tables}
DB=${DB:-/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/db/BFD.duckdb}
OUT=${OUT:-/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/summary_ecology}
PY=${PY:-/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/atlas_env/venv/bin/python}
: "${SCRATCH:?SCRATCH not set}"
export MPLCONFIGDIR="$SCRATCH/mplconfig"
mkdir -p "$OUT" "$MPLCONFIGDIR"
cd "$ROOT"
START=$(date +%s)
$PY bin/summary/summary_ecology.py --tables "$DATA" --db "$DB" --out "$OUT"
echo "elapsed_seconds=$(( $(date +%s) - START ))"
