#!/usr/bin/bash -l
#SBATCH -p short -c 2 --mem 8gb --time 1:00:00 -J kx2ft_inputs
#SBATCH -o kx2ft_inputs.%j.log
# Step 5 wrapper. Submit after step 4 finished:  sbatch --dependency=afterok:<array job id> scripts/05_build_inputs.sh
set -euo pipefail
BENCH=$(readlink -f "$HOME/projects/BFD_dashboard/bench/kx2_ku_family")
# Layer: 'annotated' (default; references + annotated candidates) or 'rescue' (adds the miniprot models): sbatch scripts/05_build_inputs.sh rescue
/usr/bin/python3.12 "$BENCH/scripts/05_build_inputs.py" "$BENCH" "${1:-annotated}"
