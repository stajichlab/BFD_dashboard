#!/usr/bin/bash -l
#SBATCH -p epyc -c 4 --mem 48G --time 1:30:00 -J atlas_J10_hs_render
#SBATCH -o /bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/atlas_data/logs/J10_hs_render.%j.log
# Render the size-transition pages (HS_G1, HS_G3, HS_P2, HS_S1, HS_C1, HS_C3) from A11 + small intermediates.
# ONLY env var selects views, e.g. ONLY=g1,g3 sbatch J10_hs_render.sh
set -euo pipefail
: "${SCRATCH:?SCRATCH is not set; this job must run under SLURM on UCR HPCC}"
ATLAS=/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard
PY="$ATLAS/atlas_env/venv/bin/python"
mkdir -p "$SCRATCH/atlas_render_tmp"
"$PY" "$ATLAS/bin/atlas/render_hs.py" ${ONLY:+--only "$ONLY"}
