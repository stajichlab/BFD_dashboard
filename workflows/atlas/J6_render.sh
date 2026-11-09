#!/usr/bin/bash -l
#SBATCH -p short -c 4 --mem 48G --time 2:00:00 -J atlas_render
#SBATCH -o /bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/atlas_data/logs/J6_render.%j.log
# Render atlas views. ONLY env var selects views (comma list), e.g. ONLY=v08,v12 sbatch J6_render.sh
set -euo pipefail
ATLAS=/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard
PY="$ATLAS/atlas_env/venv/bin/python"; [ -x "$PY" ] || PY="$ATLAS/atlas_env/.pixi/envs/default/bin/python"
if ! "$PY" -c "import statsmodels, duckdb, matplotlib" 2>/dev/null; then PY=/usr/bin/python3.12; EXTRA=--skip-mixed; echo "pixi env not ready; using $PY without mixed models"; else EXTRA=""; fi
mkdir -p "${SCRATCH:?}/atlas_render_tmp"
sed -i "s|SET temp_directory='[^']*'|SET temp_directory='${SCRATCH}/atlas_render_tmp'|" "$ATLAS/bin/atlas/render_views.py"
"$PY" "$ATLAS/bin/atlas/render_views.py" ${ONLY:+--only "$ONLY"} $EXTRA ${NOINDEX:+--no-index}
