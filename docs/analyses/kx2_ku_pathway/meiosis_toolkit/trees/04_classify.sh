#!/usr/bin/bash -l
#SBATCH -p short -c 2 --mem 8gb --time 1:00:00 -J kx2mt_classify
#SBATCH -o kx2mt_classify.%j.log
# Step 4 wrapper. Submit from the BENCH folder after the trees:  BENCH=... sbatch --dependency=afterok:<tree array id> --export=ALL scripts_trees/04_classify.sh
set -euo pipefail
: "${BENCH:?set BENCH}"; BENCH=$(readlink -f "$BENCH")
BB=$(readlink -f "$HOME/projects/BFD_dashboard/bench/kickx_tree/taxa/backbone_taxa.tsv")
/usr/bin/python3.12 "$BENCH/scripts_trees/04_classify.py" "$BENCH" "$BB" "$BENCH/scripts_trees/coemansia_strains.txt"
