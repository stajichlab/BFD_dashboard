#!/usr/bin/bash -l
#SBATCH -p short -c 2 --mem 8gb --time 1:00:00 -J kx2ft_validate
#SBATCH -o kx2ft_validate.%j.log
# Step 5v wrapper. Submit after step 4 (whatever its outcome):  sbatch --dependency=afterany:<array job id> scripts/05v_validate_models.sh
set -euo pipefail
BENCH=$(readlink -f "$HOME/projects/BFD_dashboard/bench/kx2_ku_family")
/usr/bin/python3.12 "$BENCH/scripts/05v_validate_models.py" "$BENCH"
