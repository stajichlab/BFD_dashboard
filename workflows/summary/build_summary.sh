#!/usr/bin/bash -l
#SBATCH -p epyc -c 16 --mem 64G --time 4:00:00 -J qs_build
#SBATCH -o qs_build.%j.log
# Build the small per-genome tables for the Quarto summary from the local BFD.duckdb.
# Run from a bench copy of the repo, e.g. ~/projects/BFD_dashboard/bench/quarto_summary/code
set -euo pipefail
BENCH=${BENCH:-$HOME/projects/BFD_dashboard/bench/quarto_summary}
OUT=${OUT:-$BENCH/data}
PY=${PY:-/usr/bin/python3.12}
cd "$BENCH/code"
$PY bin/summary/build_summary_tables.py "$OUT" --threads "${SLURM_CPUS_PER_TASK:-16}" --memory 48GB
ls -l "$OUT"
