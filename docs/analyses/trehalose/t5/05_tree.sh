#!/usr/bin/bash -l
#SBATCH -p epyc -c 32 --mem 96gb --time 7-00:00:00 -J t5_tree
#SBATCH -o t5_tree.%j.log
# T5 step 5: align, trim and build the GH65 maximum-likelihood tree. Submit from the bench folder after step 4:
#   BENCH=~/projects/BFD_dashboard/bench/trehalose_T5 sbatch --export=ALL code/docs/analyses/trehalose/t5/05_tree.sh
# Choices (fixed before the run, not tuned to the result):
#   mafft 7.505 --auto           the set is diverged and of uneven length (references 650-1,200 aa, fungal proteins up to about 1,300 aa)
#   trimAl 1.5.1 -gt 0.4 -resoverlap 0.5 -seqoverlap 50   keep columns present in >= 40% of sequences, drop sequences covering < 50% of the
#                                alignment; the setting used for the RecA tree after -automated1 left too few columns
#   IQ-TREE 3.0.1 -m Q.PFAM+R8   a fixed model, as for the RecA tree: ModelFinder over thousands of sequences is not practical, and Q.PFAM +
#                                a free-rate mixture is the model ModelFinder chose for the earlier diverged family trees; 1,000 ultrafast
#                                bootstrap and 1,000 SH-aLRT replicates (the support pair that fam_lib.classify reads)
# Temp files go to $SCRATCH; results are copied back to $BENCH/trees. Sequences dropped by trimAl are recorded in the log.
set -euo pipefail
: "${SCRATCH:?SCRATCH is not set; this job must run under SLURM on UCR HPCC}"
: "${BENCH:?set BENCH}"
BENCH=$(readlink -f "$BENCH"); T="$BENCH/trees"
module load mafft/7.505 trimal/1.5.1 iqtree/3.0.1
w="$SCRATCH/t5_tree"; mkdir -p "$w"; cd "$w"
in="$T/gh65_input.faa"; test -s "$in" || { echo "missing $in" >&2; exit 2; }
{ echo "mafft: $(mafft --version 2>&1 | head -1)"; echo "trimal: $(trimal --version 2>&1 | head -2 | tail -1)"; echo "iqtree: $(iqtree3 --version 2>&1 | head -1)"; } > "$T/gh65_versions.txt"
echo "sequences: $(grep -c '^>' "$in")"
mafft --auto --thread 32 "$in" > aln.fa 2> mafft.log
trimal -in aln.fa -out trim.fa -gt 0.4 -resoverlap 0.5 -seqoverlap 50 -fasta
echo "alignment columns: $(awk '/^>/{n++; next} n==1{c+=length($0)} END{print c}' aln.fa) -> trimmed: $(awk '/^>/{n++; next} n==1{c+=length($0)} END{print c}' trim.fa)"
echo "sequences after trimming: $(grep -c '^>' trim.fa)"
cp aln.fa "$T/gh65.aln.fa"; cp trim.fa "$T/gh65.trim.fa"
iqtree3 -s trim.fa -m Q.PFAM+R8 -B 1000 --alrt 1000 -T 32 --prefix gh65 -seed 20261113 --mem 90G
cp gh65.treefile gh65.iqtree gh65.log "$T/"
echo done
