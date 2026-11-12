#!/usr/bin/bash -l
#SBATCH -p short -c 8 --mem 32gb --time 2:00:00 -J kx2ft_cand
#SBATCH -o kx2ft_cand.%j.log
# Step 2 job wrapper. Submit from the BENCH folder:  sbatch $BENCH/scripts/02_candidates.sh
# Inputs: the backbone taxa list written by the Kickxellomycota tree job (read-only) and strain_genomes.tsv, a copy of
# ../tblastn/genomes.tsv (the 31 worklist genomes, which include the extra Coemansia strains).
set -euo pipefail
BENCH=$(readlink -f "$HOME/projects/BFD_dashboard/bench/kx2_ku_family")
BB=$(readlink -f "$HOME/projects/BFD_dashboard/bench/kickx_tree/taxa/backbone_taxa.tsv")
/usr/bin/python3.12 "$BENCH/scripts/02_candidates.py" "$BENCH" "$BB" "$BENCH/scripts/strain_genomes.tsv"
