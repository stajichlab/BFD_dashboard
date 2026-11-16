#!/usr/bin/bash -l
#SBATCH -p epyc -c 16 --mem 128G --time 12:00:00 -J atlas_J12
#SBATCH -o kx2_j12.%j.log
# J12: end-to-end rebuild of Analysis E (H1-H4) from the local BFD DuckDB into an empty scratch database.
# Reads BFD.duckdb (read-only, shared). Writes only under $BENCH (a purgeable folder in the user's project, not the shared group folder):
#   $BENCH/out_<jobid>/{run_log.tsv, results/*.tsv, table_counts.tsv, table_fingerprint.tsv}
# The scratch DuckDB file and DuckDB temp files go to node-local $SCRATCH and disappear when the job ends.
#
# Submit from $BENCH after copying the code there (so the run is pinned to a known copy of the repo files):
#   BENCH=~/projects/BFD_dashboard/bench/kx2_j12 sbatch --export=ALL workflows/atlas/J12_kx2_rebuild.sh
# Optional: BFD_DB=/path/BFD.duckdb  J12_STAGES=upstream,h1,h2,h3,h4  J12_CODE=<dir holding bin/ lib/ docs/>
# Test: /usr/bin/python3.12 -m pytest bin/atlas/test_j12_kx2_rebuild.py   (small fixture, no HPCC data)
set -euo pipefail
: "${SCRATCH:?SCRATCH is not set; this job must run under SLURM on UCR HPCC}"
: "${BENCH:?set BENCH to the output folder, e.g. ~/projects/BFD_dashboard/bench/kx2_j12}"
BENCH=$(readlink -f "$BENCH")
CODE=$(readlink -f "${J12_CODE:-$BENCH/code}")
BFD_DB=${BFD_DB:-/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/db/BFD.duckdb}
OUT="$BENCH/out_${SLURM_JOB_ID:-manual}"
mkdir -p "$OUT" "$SCRATCH/duckdb_tmp"
rm -f "$SCRATCH/j12.duckdb"          # always start from an empty scratch database
/usr/bin/python3.12 "$CODE/bin/atlas/j12_kx2_rebuild.py" \
  --sql-root "$CODE/docs/analyses/kx2_kickx_reduction" \
  --lib "$CODE/lib" \
  --bfd-db "$BFD_DB" \
  --scratch-db "$SCRATCH/j12.duckdb" \
  --tmp "$SCRATCH/duckdb_tmp" \
  --threads 16 --memory 110GB \
  --stages "${J12_STAGES:-upstream,h1,h2,h3,h4}" \
  --out "$OUT"
