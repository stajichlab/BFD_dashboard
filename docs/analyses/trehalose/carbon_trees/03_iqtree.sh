#!/usr/bin/bash -l
#SBATCH -p epyc -c 32 --mem 96gb --time 5-00:00:00 -J ct_iq
#SBATCH --array=1-4
#SBATCH -o ct_iq.%A_%a.log
# IQ-TREE maximum-likelihood trees for the four families of at most 2,300 sequences, after 02_trees.sh. Submit from the bench folder:
#   BENCH=$PWD sbatch code/docs/analyses/trehalose/carbon_trees/03_iqtree.sh
# IQ-TREE 3.0.1 -m Q.PFAM+R8, 1,000 ultrafast bootstrap + 1,000 SH-aLRT, fixed seed (as T5 and as fixed in README.md).
set -euo pipefail
: "${SCRATCH:?SCRATCH is not set; this job must run under SLURM on UCR HPCC}"
: "${BENCH:?set BENCH}"
BENCH=$(readlink -f "$BENCH"); T="$BENCH/trees"
FAMS=(GSY GPH GDB GPD)
fam=${FAMS[$((SLURM_ARRAY_TASK_ID-1))]}
module load iqtree/3.0.1
w="$SCRATCH/ctiq_$fam"; mkdir -p "$w"; cd "$w"
test -s "$T/$fam.trim.fa" || { echo "missing $T/$fam.trim.fa (run 02_trees.sh first)" >&2; exit 2; }
cp "$T/$fam.trim.fa" trim.fa
echo "$fam: $(grep -c '^>' trim.fa) sequences, $(awk '/^>/{n++; next} n==1{c+=length($0)} END{print c}' trim.fa) columns"
iqtree3 -s trim.fa -m Q.PFAM+R8 -B 1000 --alrt 1000 -T 32 --prefix "$fam" -seed 20261116 --mem 90G
cp "$fam.treefile" "$fam.iqtree" "$fam.log" "$T/"
echo done
