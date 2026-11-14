#!/usr/bin/bash -l
#SBATCH -p epyc -c 16 --mem 100G --time 3:00:00 -J atlas_J13_trehalose
#SBATCH -o trehalose_J13.%j.log
# Atlas job J13: rebuild the trehalose survey (analysis G, docs/analyses/trehalose/) from the BFD parquet tables.
#   step 0  fixture tests (BH q-values, class-stratified co-absence null)
#   step 1  th0_inputs.sql  -> $TH_BENCH/data/th_*.csv      (DuckDB; restartable: existing outputs are skipped)
#   step 2  th1..th5        -> $TH_BENCH/data/th_*.tsv, run_log_*.txt, $TH_BENCH/figs/*.png
# Submit from the bench folder, with a checkout of the repo that contains this file:
#   mkdir -p ~/projects/BFD_dashboard/bench/trehalose_J13 && cd $_
#   TH_REPO=/path/to/checkout sbatch --export=ALL $TH_REPO/workflows/atlas/J13_trehalose.sh
# Outputs go to a purgeable bench folder, NOT the shared atlas_data directory (the shared folder is read-only here:
# it is read for A7/pfam_genome_matrix.parquet and the Fungi_BFD_runs/tables parquet files).
# Set TH_FORCE=1 to rebuild the inputs.
set -euo pipefail
: "${SCRATCH:?SCRATCH is not set; this job must run under SLURM on UCR HPCC}"
: "${TH_REPO:?set TH_REPO to a checkout containing this branch}"
TH_BENCH=${TH_BENCH:-$HOME/projects/BFD_dashboard/bench/trehalose_J13}
TH_BENCH=$(readlink -f "$TH_BENCH")
export TH_DATA="$TH_BENCH/data" TH_FIGS="$TH_BENCH/figs"
mkdir -p "$TH_DATA" "$TH_FIGS"
PY=/usr/bin/python3.12
S="$TH_REPO/docs/analyses/trehalose/scripts"

echo "[J13] repo $(git -C "$TH_REPO" rev-parse --short HEAD 2>/dev/null || echo unknown), bench $TH_BENCH"
cd "$S"
$PY test_stats_util.py
$PY "$TH_REPO/bin/atlas/run_duckdb_job.py" --sql "$S/th0_inputs.sql" --atlas "$TH_BENCH" \
    --tmp "$SCRATCH/duckdb_tmp" --threads 16 --memory 90GB ${TH_FORCE:+--force}
for s in th1_survey th2_stats th3_coloss th4_predict th5_qc; do
  echo "[J13] $s"
  $PY "$S/$s.py" 2>&1 | tee "$TH_DATA/run_log_${s%%_*}_${s#*_}.txt"
done
echo "[J13] done"
