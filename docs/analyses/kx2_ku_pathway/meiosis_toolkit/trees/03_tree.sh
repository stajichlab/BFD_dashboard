#!/usr/bin/bash -l
#SBATCH -p epyc -c 16 --mem 32gb --time 24:00:00 -J kx2mt_tree
#SBATCH --array=1-4
#SBATCH -o kx2mt_tree.%A_%a.log
# Step 3: align, trim and build a maximum-likelihood tree for one family (array task 1 = muts, 2 = mutl, 3 = kleisin, 4 = horma).
# Submit from the BENCH folder after steps 1-2:  BENCH=... sbatch --export=ALL 03_tree.sh
#   mafft 7.505 --auto            members are long, diverged and of uneven length; --auto chooses an accurate strategy for the set size
#   trimAl 1.5.1 -gt 0.4 -resoverlap 0.5 -seqoverlap 50   the GENTLE mode that worked for the RecA re-trim in ../../family_trees/ (its
#                                 first tree used -automated1 and kept only 41 columns, which made it unusable); the alignment width is
#                                 checked below and the job STOPS if fewer than 150 columns remain, so a thin tree is never built
#   IQ-TREE 3.0.1 -m LG+G4        one fixed, well-tested model for speed (the RecA tree with ModelFinder took 1.5 h for 362 sequences);
#                                 ultrafast bootstrap (1000) and SH-aLRT (1000), the support convention that fam_lib.classify reads
set -euo pipefail
: "${SCRATCH:?SCRATCH is not set; this job must run under SLURM on UCR HPCC}"; : "${BENCH:?set BENCH}"
BENCH=$(readlink -f "$BENCH")
fam=$(echo "muts mutl kleisin horma" | cut -d' ' -f"${SLURM_ARRAY_TASK_ID}")
module load mafft/7.505 trimal/1.5.1 iqtree/3.0.1
w="$SCRATCH/kx2mt_$fam"; mkdir -p "$w"; cd "$w"
in="$BENCH/trees/${fam}_input.faa"; test -s "$in" || { echo "missing $in" >&2; exit 2; }
{ echo "mafft: $(mafft --version 2>&1 | head -1)"; echo "trimal: $(trimal --version 2>&1 | head -2 | tail -1)"; echo "iqtree: $(iqtree3 --version 2>&1 | head -1)"; } > "$BENCH/trees/${fam}_versions.txt"
echo "$fam: $(grep -c '^>' "$in") sequences"
mafft --auto --thread 16 "$in" > aln.fa 2> mafft.log
trimal -in aln.fa -out trim.fa -gt 0.4 -resoverlap 0.5 -seqoverlap 50 -fasta
cols=$(awk '/^>/{n++; next} n==1{c+=length($0)} END{print c}' trim.fa)
echo "alignment columns: $(awk '/^>/{n++; next} n==1{c+=length($0)} END{print c}' aln.fa) -> trimmed: $cols; sequences kept: $(grep -c '^>' trim.fa)"
if [ "$cols" -lt 150 ]; then echo "STOP: only $cols columns after trimming" >&2; exit 7; fi
iqtree3 -s trim.fa -m LG+G4 -B 1000 --alrt 1000 -T 16 --prefix "$fam" -seed 20261112
cp aln.fa "$BENCH/trees/${fam}.aln.fa"; cp trim.fa "$BENCH/trees/${fam}.trim.fa"
cp "$fam.treefile" "$fam.iqtree" "$fam.log" "$BENCH/trees/"
echo done
