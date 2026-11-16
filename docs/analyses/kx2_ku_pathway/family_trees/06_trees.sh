#!/usr/bin/bash -l
#SBATCH -p epyc -c 16 --mem 32gb --time 24:00:00 -J kx2ft_tree
#SBATCH --array=1-3
#SBATCH -o kx2ft_tree.%A_%a.log
# Step 6: align, trim and build a maximum-likelihood tree for one family (array task 1 = lig, 2 = reca, 3 = polx).
# Submit after step 5:  sbatch scripts/06_trees.sh
#   mafft 7.505 --auto           family members are long, diverged and of uneven length; --auto picks an accurate strategy for the
#                                size of the set (the sets are a few hundred sequences)
#   trimAl 1.5.1 -automated1     removes poorly aligned columns using a heuristic tuned for maximum likelihood trees; the long variable
#                                regions of these proteins (ligase N- and C-terminal regions, Pol X N-terminus) would otherwise dominate
#   IQ-TREE 3.0.1 -m MFP         ModelFinder over LG, WAG, JTT and Q.pfam with rate-heterogeneity variants (BIC), then the best model;
#                                ultrafast bootstrap (1000) and SH-aLRT (1000), the support convention that 07_classify.py reads
# Work is done in $SCRATCH; results are copied back to $BENCH/trees.
set -euo pipefail
: "${SCRATCH:?SCRATCH is not set; this job must run under SLURM on UCR HPCC}"
BENCH=$(readlink -f "$HOME/projects/BFD_dashboard/bench/kx2_ku_family")
fam=$(echo "lig reca polx" | cut -d' ' -f"${SLURM_ARRAY_TASK_ID}")
LAYER=${LAYER:-annotated}                        # annotated | rescue (set with: LAYER=rescue sbatch --export=ALL scripts/06_trees.sh)
T="$BENCH/trees_$LAYER"
module load mafft/7.505 trimal/1.5.1 iqtree/3.0.1
w="$SCRATCH/kx2ft_${LAYER}_$fam"; mkdir -p "$w"; cd "$w"
in="$T/${fam}_input.faa"; test -s "$in" || { echo "missing $in" >&2; exit 2; }
{ echo "mafft: $(mafft --version 2>&1 | head -1)"; echo "trimal: $(trimal --version 2>&1 | head -2 | tail -1)"; echo "iqtree: $(iqtree3 --version 2>&1 | head -1)"; } > "$T/${fam}_versions.txt"
echo "$fam: $(grep -c '^>' "$in") sequences"
mafft --auto --thread 16 "$in" > aln.fa 2> mafft.log
trimal -in aln.fa -out trim.fa -automated1 -fasta
echo "alignment columns: $(awk '/^>/{n++; next} n==1{c+=length($0)} END{print c}' aln.fa) -> trimmed: $(awk '/^>/{n++; next} n==1{c+=length($0)} END{print c}' trim.fa)"
iqtree3 -s trim.fa -m MFP -mset LG,WAG,JTT,Q.pfam -B 1000 --alrt 1000 -T 16 --prefix "$fam" -seed 20261111
cp aln.fa "$T/${fam}.aln.fa"; cp trim.fa "$T/${fam}.trim.fa"
cp "$fam.treefile" "$fam.iqtree" "$fam.log" "$T/"
echo done
