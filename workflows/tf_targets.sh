#!/usr/bin/bash -l
#SBATCH -p epyc -c 8 --mem 40G --time 4:00:00 -J tf_targets -o %x.%j.log
# DAP-seq target sets (Lax et al. 2025, R. microsporus ATCC 11559): JGI-to-BFD protein map, target conservation, profiles.
# Submit from the repo root: sbatch workflows/tf_targets.sh
# Steps (each skips when its output exists): prep, map (RBH 90/0.8 and 70/0.7), repdb (genus-representative MMseqs2 DB,
# about 16 million proteins), conserve (all FBF4BF29 proteins vs that DB), tables, figures.
# Large files stay in the staging dir. MMseqs2 temp files go to node-local $SCRATCH. Never BASH_SOURCE here.
set -euo pipefail
ROOT=${ROOT:-${SLURM_SUBMIT_DIR:?submit from the repo root or set ROOT}}
STAGE=${STAGE:-/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/tf_targets}
SUMMARY=${SUMMARY:-/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/summary_tables}
PY=${PY:-/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/atlas_env/venv/bin/python}
: "${SCRATCH:?}"
module load MMseqs2/17-b804f
mkdir -p "$STAGE"
cd "$ROOT"
export SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-8}
T=bin/tf/targets.py
for step in prep map repdb conserve tables figures; do
  echo "== step $step $(date +%T)"
  $PY $T "$step" --stage "$STAGE" --scratch "$SCRATCH" --summary-dir "$SUMMARY"
done
ls -l "$STAGE/tables_tf/" "$STAGE/figures/tf/" "$STAGE/genus_reps_db.info.tsv"
cat "$STAGE/genus_reps_db.info.tsv"
