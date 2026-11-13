#!/usr/bin/bash -l
#SBATCH -p short -c 8 --mem 32gb --time 2:00:00 -J kx2mt_inputs
#SBATCH -o kx2mt_inputs.%j.log
# Step 2 wrapper. Submit from the BENCH folder after step 1:  BENCH=... sbatch --export=ALL scripts_trees/02_build_inputs.sh
set -euo pipefail
: "${BENCH:?set BENCH}"; BENCH=$(readlink -f "$BENCH")
BB=$(readlink -f "$HOME/projects/BFD_dashboard/bench/kickx_tree/taxa/backbone_taxa.tsv")
/usr/bin/python3.12 "$BENCH/scripts_trees/02_build_inputs.py" "$BENCH" "$BB" "$BENCH/scripts_trees/coemansia_strains.txt"
