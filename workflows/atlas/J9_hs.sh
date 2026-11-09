#!/usr/bin/bash -l
#SBATCH -p epyc -c 16 --mem 120G --time 3:00:00 -J atlas_J9_hs
#SBATCH -o /bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/atlas_data/logs/J9_hs.%j.log
# Size-transition views: build A11 intermediates (DuckDB). Restartable: existing outputs are skipped (run_duckdb_job.py).
set -euo pipefail
: "${SCRATCH:?SCRATCH is not set; this job must run under SLURM on UCR HPCC}"
ATLAS=/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard
/usr/bin/python3.12 "$ATLAS/bin/atlas/run_duckdb_job.py" \
  --sql "$ATLAS/bin/atlas/sql/J9_hs.sql" \
  --atlas "$ATLAS/atlas_data" \
  --tmp "$SCRATCH/duckdb_tmp" \
  --threads 16 --memory 100GB ${J9_FORCE:+--force}
