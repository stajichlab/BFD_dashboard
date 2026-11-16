#!/bin/bash
# Submit analysis C on the UCR HPCC as four chained SLURM jobs: prep -> (rf, phylo in parallel) -> report.
# Usage (from a checkout of the repo on the HPCC):
#   AO_BASE=~/projects/BFD_dashboard/bench/antarctolithica_C/work docs/analyses/antarctolithica/run_hpcc.sh
# Optional: BFD_DATA (shared BFD root, read only), AO_PHYLO_TAXA (taxa.tsv fixing the tree taxa; default the committed results/taxa.tsv).
# One-time setup: (cd docs/analyses/antarctolithica && pixi install)
set -euo pipefail
export AO_REPO=$(cd "$(dirname "$0")/../../.." && pwd)
: "${AO_BASE:?set AO_BASE to the working directory}"
mkdir -p "$AO_BASE/logs" "$AO_BASE/results"
export AO_BASE=$(cd "$AO_BASE" && pwd) BFD_DATA=${BFD_DATA:-/bigdata/stajichlab/shared/projects/BFD}
export AO_COMMIT=$(git -C "$AO_REPO" rev-parse --short HEAD 2>/dev/null || echo "${AO_COMMIT:-unknown}")
SB=$AO_REPO/docs/analyses/antarctolithica/run_hpcc.sbatch
cd "$AO_BASE"
p=$(sbatch --parsable --export=ALL --job-name=aoC_prep  --cpus-per-task=16 --mem=48G  --time=8:00:00  "$SB" prep)
r=$(sbatch --parsable --export=ALL --job-name=aoC_rf    --cpus-per-task=16 --mem=32G  --time=24:00:00 --dependency=afterok:$p "$SB" rf)
f=$(sbatch --parsable --export=ALL --job-name=aoC_phylo --cpus-per-task=16 --mem=16G  --time=24:00:00 --dependency=afterok:$p "$SB" phylo)
o=$(sbatch --parsable --export=ALL --job-name=aoC_report --cpus-per-task=4 --mem=16G  --time=2:00:00  --dependency=afterok:$r:$f "$SB" report)
echo "commit $AO_COMMIT  AO_BASE $AO_BASE"
echo "jobs: prep=$p rf=$r phylo=$f report=$o"
