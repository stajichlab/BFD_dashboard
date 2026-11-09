#!/usr/bin/bash -l
#SBATCH -p epyc -c 16 --mem 128G --time 6:00:00 -J atlas_J3
#SBATCH -o /bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/atlas_data/logs/J3.%j.log
# Atlas intermediate job J3. Restartable: existing outputs are skipped (see run_duckdb_job.py).
set -euo pipefail
: "${SCRATCH:?SCRATCH is not set; this job must run under SLURM on UCR HPCC}"
ATLAS=/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard
/usr/bin/python3.12 "$ATLAS/bin/atlas/run_duckdb_job.py" \
  --sql "$ATLAS/bin/atlas/sql/J3_protein_core.sql" \
  --atlas "$ATLAS/atlas_data" \
  --tmp "$SCRATCH/duckdb_tmp" \
  --threads 16 --memory 112GB
