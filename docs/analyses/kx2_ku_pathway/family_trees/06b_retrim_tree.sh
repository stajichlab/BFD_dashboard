#!/usr/bin/bash -l
#SBATCH -p epyc -c 16 --mem 32gb --time 24:00:00 -J kx2ft_retrim
#SBATCH -o kx2ft_retrim.%j.log
# Step 6b (a deviation from step 6, made once): re-trim and re-run the tree of one family on one layer with a gentler trimAl setting.
# Why: trimAl -automated1 reduced the RecA/Rad51 alignment (362 sequences) to 41 columns, which cannot separate Rad51, Dmc1, Rad55 and
# Rad57. Here columns present in at least 40% of sequences are kept (-gt 0.4), after dropping sequences that are mostly gaps (-resoverlap 0.5
# -seqoverlap 50 removes sequences that overlap the rest of the alignment by less than half of the columns).
# Usage:  LAYER=annotated FAM=reca sbatch --export=ALL scripts/06b_retrim_tree.sh      (the old tree is kept as <fam>.treefile.automated1)
set -euo pipefail
: "${SCRATCH:?}"; : "${FAM:?set FAM}"
BENCH=$(readlink -f "$HOME/projects/BFD_dashboard/bench/kx2_ku_family")
LAYER=${LAYER:-annotated}; T="$BENCH/trees_$LAYER"; fam=$FAM
module load trimal/1.5.1 iqtree/3.0.1
w="$SCRATCH/kx2ft_retrim_${LAYER}_$fam"; mkdir -p "$w"; cd "$w"
cp "$T/$fam.aln.fa" aln.fa
trimal -in aln.fa -out trim.fa -gt 0.4 -resoverlap 0.5 -seqoverlap 50 -fasta
echo "$fam $LAYER: sequences $(grep -c '^>' aln.fa) -> $(grep -c '^>' trim.fa) after overlap filter; columns $(awk '/^>/{n++; next} n==1{c+=length($0)} END{print c}' aln.fa) -> $(awk '/^>/{n++; next} n==1{c+=length($0)} END{print c}' trim.fa)"
for ext in treefile iqtree trim.fa; do [ -e "$T/$fam.$ext" ] && cp "$T/$fam.$ext" "$T/$fam.$ext.automated1"; done
iqtree3 -s trim.fa -m MFP -mset LG,WAG,JTT,Q.pfam -B 1000 --alrt 1000 -T 16 --prefix "$fam" -seed 20261111
cp trim.fa "$T/$fam.trim.fa"; cp "$fam.treefile" "$fam.iqtree" "$fam.log" "$T/"
echo done
