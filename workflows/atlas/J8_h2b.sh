#!/usr/bin/bash -l
#SBATCH -p epyc -c 16 --mem 96G --time 3:00:00 -J atlas_J8
#SBATCH -o /bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/atlas_data/logs/J8.%j.log
# H2b: build A10 intermediates (DuckDB) then render docs/summary_atlas/H2b_compaction_extremes.md. Restartable (skips existing A10 outputs).
set -euo pipefail
: "${SCRATCH:?SCRATCH is not set}"
ATLAS=/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard
PY="$ATLAS/atlas_env/venv/bin/python"
if [ "${H2B_SKIP_BUILD:-0}" != "1" ]; then
  /usr/bin/python3.12 "$ATLAS/bin/atlas/h2b_build.py" --tmp "$SCRATCH/duckdb_tmp" --threads 16 --memory 80GB ${H2B_FORCE:+--force}
fi
mkdir -p "$SCRATCH/atlas_render_tmp"
"$PY" "$ATLAS/bin/atlas/render_h2b.py"
