#!/usr/bin/bash -l
#SBATCH -p epyc -c 16 --mem 64gb --time 1-00:00:00 -J ct_fast
#SBATCH --array=1-8
#SBATCH -o ct_fast.%A_%a.log
# First-pass trees for the carbon-storage families, one array task per family. Submit from the bench folder after 01_extract_domains.py:
#   BENCH=$PWD sbatch code/docs/analyses/trehalose/carbon_trees/02_trees.sh
# mafft 7.505 --auto; trimAl 1.5.1 -gt 0.4 -resoverlap 0.5 -seqoverlap 50; FastTree 2.1.11 -lg -gamma (as fixed in README.md).
# Temp files go to $SCRATCH; alignments, trimmed alignments and trees are copied to $BENCH/trees.
set -euo pipefail
: "${SCRATCH:?SCRATCH is not set; this job must run under SLURM on UCR HPCC}"
: "${BENCH:?set BENCH}"
BENCH=$(readlink -f "$BENCH"); T="$BENCH/trees"; mkdir -p "$T"
FAMS=(TPS TPP NTH UGP GSY GPH GDB GPD)
fam=${FAMS[$((SLURM_ARRAY_TASK_ID-1))]}
module load mafft/7.505 trimal/1.5.1 fasttree/2.1.11
in="$BENCH/fam/$fam.faa"; test -s "$in" || { echo "missing $in" >&2; exit 2; }
w="$SCRATCH/ct_$fam"; mkdir -p "$w"; cd "$w"
echo "$fam: $(grep -c '^>' "$in") sequences"
mafft --auto --thread 16 "$in" > aln.fa 2> mafft.log
trimal -in aln.fa -out trim.fa -gt 0.4 -resoverlap 0.5 -seqoverlap 50 -fasta
echo "alignment columns: $(awk '/^>/{n++; next} n==1{c+=length($0)} END{print c}' aln.fa) -> trimmed: $(awk '/^>/{n++; next} n==1{c+=length($0)} END{print c}' trim.fa)"
echo "sequences after trimming: $(grep -c '^>' trim.fa)"
cp trim.fa "$T/$fam.trim.fa"; cp aln.fa "$T/$fam.aln.fa"
FastTree -lg -gamma < trim.fa > "$T/$fam.fasttree.nwk" 2> "$T/$fam.fasttree.log"
{ echo "mafft: $(mafft --version 2>&1 | head -1)"; echo "trimal: $(trimal --version 2>&1 | head -2 | tail -1)"; } > "$T/$fam.versions.txt"
echo done
