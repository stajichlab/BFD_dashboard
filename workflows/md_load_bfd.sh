#!/usr/bin/bash -l
#SBATCH -p epyc -c 16 --mem 48G --time 24:00:00 -J md_load_bfd
#SBATCH -o md_load_bfd.%j.log
# Load BFD parquet tables into MotherDuck (database "bfd" by default).
# Restartable: resubmit after a failure/timeout and finished tables are skipped.
#
# One-time setup (on HPCC, not in this script):
#   mkdir -p ~/.config/motherduck && chmod 700 ~/.config/motherduck
#   cat > ~/.config/motherduck/token      # paste token from MotherDuck Settings > Access Tokens, then Ctrl-D
#   chmod 600 ~/.config/motherduck/token
#
# Submit:   sbatch workflows/md_load_bfd.sh            (all tables except swissprot_hits)
#           sbatch workflows/md_load_bfd.sh --only species asm_stats pfam
set -euo pipefail
: "${SCRATCH:?SCRATCH is not set; run under SLURM on UCR HPCC}"
REPO=/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard
TABLES=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables

# Fail fast if this node cannot reach MotherDuck.
curl -s -m 20 -o /dev/null https://api.motherduck.com/ \
  || { echo "No outbound HTTPS to api.motherduck.com from $(hostname)"; exit 2; }

/usr/bin/python3.12 "$REPO/bin/md_load_bfd.py" \
  --tables "$TABLES" --db bfd \
  --threads "${SLURM_CPUS_PER_TASK:-16}" --memory 40GB \
  --tmp "$SCRATCH/md_duckdb_tmp" "$@"
