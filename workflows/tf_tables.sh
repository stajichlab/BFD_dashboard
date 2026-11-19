#!/usr/bin/bash -l
#SBATCH -p epyc -c 4 --mem 32G --time 2:00:00 -J tf_tables -o %x.%j.log
# Build the per-genome TF tables (tf_family_genome, tf_summary_genome, tf_paper_family_genome) for ALL genomes.
# Submit from the repo root: sbatch workflows/tf_tables.sh
# Reads Fungi_BFD_runs/tables/{tf_inventory,pfam,gene_proteins}.parquet and the summary tables.
# DuckDB temp files go to node-local $SCRATCH. Outputs stay on /bigdata.
set -euo pipefail
ROOT=${ROOT:-${SLURM_SUBMIT_DIR:?submit from the repo root or set ROOT}}
TABLES=${TABLES:-/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables}
OUT=${OUT:-/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/summary_tables}
PY=${PY:-/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/atlas_env/venv/bin/python}
: "${SCRATCH:?}"
mkdir -p "$OUT"
cd "$ROOT"
$PY bin/summary/build_tf_tables.py "$TABLES" "$OUT" "$OUT" --threads "${SLURM_CPUS_PER_TASK:-4}" --memory 28GB --tmp "$SCRATCH"
ls -l "$OUT/tf_family_genome.parquet" "$OUT/tf_summary_genome.parquet" "$OUT/tf_paper_family_genome.parquet"
