#!/usr/bin/bash -l
#SBATCH -p epyc -c 16 --mem 64gb --time 48:00:00 -J kx2ra_tree
#SBATCH -o kx2ra_tree.%j.log
# RecA/Rad51 tree for all Coemansia strains (annotated layer). Repeats the final RecA tree of ../06b_retrim_tree.sh on the larger input:
#   mafft --auto, trimAl -gt 0.4 -resoverlap 0.5 -seqoverlap 50 (columns present in >= 40% of sequences; the first -automated1 trim left 41 columns),
#   IQ-TREE 3.0.1, 1000 ultrafast bootstraps + 1000 SH-aLRT. DEVIATION: the model is fixed to Q.PFAM+R8, the ModelFinder choice for the
#   362-sequence tree (../results/annotated/reca.iqtree), instead of re-running ModelFinder on ~800 sequences; it keeps the two trees comparable.
# Submit after 05_build_inputs.sh:  sbatch --dependency=afterok:<id> scripts/06_tree_reca_all.sh
set -euo pipefail
: "${SCRATCH:?SCRATCH is not set; this job must run under SLURM on UCR HPCC}"
BENCH=$(readlink -f "$HOME/projects/BFD_dashboard/bench/kx2_reca_all")
T="$BENCH/trees_annotated"; fam=reca
module load mafft/7.505 trimal/1.5.1 iqtree/3.0.1
w="$SCRATCH/kx2ra_tree"; mkdir -p "$w"; cd "$w"
in="$T/${fam}_input.faa"; test -s "$in" || { echo "missing $in" >&2; exit 2; }
{ echo "mafft: $(mafft --version 2>&1 | head -1)"; echo "trimal: $(trimal --version 2>&1 | head -2 | tail -1)"; echo "iqtree: $(iqtree3 --version 2>&1 | head -1)"; } > "$T/${fam}_versions.txt"
echo "$fam: $(grep -c '^>' "$in") sequences"
mafft --auto --thread 16 "$in" > aln.fa 2> mafft.log
trimal -in aln.fa -out trim.fa -gt 0.4 -resoverlap 0.5 -seqoverlap 50 -fasta
echo "sequences $(grep -c '^>' aln.fa) -> $(grep -c '^>' trim.fa); columns $(awk '/^>/{n++; next} n==1{c+=length($0)} END{print c}' aln.fa) -> $(awk '/^>/{n++; next} n==1{c+=length($0)} END{print c}' trim.fa)"
iqtree3 -s trim.fa -m Q.PFAM+R8 -B 1000 --alrt 1000 -T 16 --prefix "$fam" -seed 20261111
cp aln.fa "$T/${fam}.aln.fa"; cp trim.fa "$T/${fam}.trim.fa"
cp "$fam.treefile" "$fam.iqtree" "$fam.log" "$T/"
echo done
