#!/usr/bin/bash -l
#SBATCH -p epyc -c 8 --mem 32G --time 3:00:00 -J tf_ortho -o %x.%j.log
# Ortholog tracking for WC-2, Zap1-like and bZip5 among the genus-representative proteomes.
# Submit from the repo root: sbatch workflows/tf_orthologs.sh
# Steps (each skips when its output exists): seeds, candidates, search, trees, report.
# Reuses the Task 5b files in $TSTAGE (genus_reps.tsv, genus_reps_db, FBF4BF29.faa); it does not rebuild them.
# Outputs go to $STAGE on /bigdata (tables/tf, figures/tf, tree and hit files). Temp files go to node-local $SCRATCH.
# Paths come from ROOT and SLURM_SUBMIT_DIR.
set -euo pipefail
ROOT=${ROOT:-${SLURM_SUBMIT_DIR:?submit from the repo root or set ROOT}}
STAGE=${STAGE:-/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/tf_orthologs}
TSTAGE=${TSTAGE:-/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/tf_targets}
SUMMARY=${SUMMARY:-/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/summary_tables}
PY=${PY:-/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/atlas_env/venv/bin/python}
: "${SCRATCH:?}"
module load MMseqs2/17-b804f
module load mafft/7.505
module load fasttree/2.1.11
mkdir -p "$STAGE"
cd "$ROOT"
export SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-8}
COMMON=(--stage "$STAGE" --scratch "$SCRATCH")
for step in seeds candidates search; do
  echo "== step $step $(date +%T)"
  $PY bin/tf/ortho_search.py "$step" "${COMMON[@]}" --targets-stage "$TSTAGE"
done
echo "== step trees $(date +%T)"
$PY bin/tf/ortho_tree.py trees "${COMMON[@]}"
echo "== step report $(date +%T)"
$PY bin/tf/ortho_report.py report --stage "$STAGE" --summary-dir "$SUMMARY"
cat "$STAGE/counts.tsv"
ls -l "$STAGE/tables/tf/" "$STAGE/figures/tf/"
