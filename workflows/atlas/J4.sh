#!/usr/bin/bash -l
#SBATCH -p epyc -c 16 --mem 96G --time 4:00:00 -J atlas_J4
#SBATCH -o /bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/atlas_data/logs/J4.%j.log
# Atlas intermediate job J4. Restartable: existing outputs are skipped (see run_duckdb_job.py).
set -euo pipefail
: "${SCRATCH:?SCRATCH is not set; this job must run under SLURM on UCR HPCC}"
ATLAS=/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard
/usr/bin/python3.12 "$ATLAS/bin/atlas/run_duckdb_job.py" \
  --sql "$ATLAS/bin/atlas/sql/J4_function_counts.sql" \
  --atlas "$ATLAS/atlas_data" \
  --tmp "$SCRATCH/duckdb_tmp" \
  --threads 16 --memory 84GB
