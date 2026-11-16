#!/usr/bin/bash -l
#SBATCH -p short -c 4 --mem 16gb --time 1:00:00 -J kx2ft_classify
#SBATCH -o kx2ft_classify.%j.log
# Step 7 wrapper. Submit after step 6:  sbatch scripts/07_classify.sh
# scripts/tblastn_round2_calls.tsv is a copy of ../tblastn/results/round2_calls.tsv (for the Lig4 comparison).
set -euo pipefail
BENCH=$(readlink -f "$HOME/projects/BFD_dashboard/bench/kx2_ku_family")
BB=$(readlink -f "$HOME/projects/BFD_dashboard/bench/kickx_tree/taxa/backbone_taxa.tsv")
# Layer: annotated (default) or rescue:  sbatch scripts/07_classify.sh rescue
export LAYER="${1:-annotated}"
/usr/bin/python3.12 "$BENCH/scripts/07a_ligase_architecture.py" "$BENCH"
/usr/bin/python3.12 "$BENCH/scripts/07_classify.py" "$BENCH" "$BB" "$BENCH/scripts/strain_genomes.tsv" "$BENCH/scripts/tblastn_round2_calls.tsv"
