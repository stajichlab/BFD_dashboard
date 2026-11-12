#!/usr/bin/bash -l
#SBATCH -p short -c 2 --mem 8gb --time 1:00:00 -J kx2ft_extract
#SBATCH -o kx2ft_extract.%j.log
# Step 3 wrapper. Submit from the BENCH folder:  sbatch scripts/03_extract_proteins.sh
set -euo pipefail
BENCH=$(readlink -f "$HOME/projects/BFD_dashboard/bench/kx2_ku_family")
/usr/bin/python3.12 "$BENCH/scripts/03_extract_proteins.py" "$BENCH"
