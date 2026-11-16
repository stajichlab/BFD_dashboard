#!/usr/bin/bash -l
#SBATCH -p short -c 2 --mem 8gb --time 0:30:00 -J kx2ft_Jtest
#SBATCH -o kx2ft_Jtest.%j.log
# Diagnostics for miniprot 0.2-r116 (not part of the pipeline). Two tests:
#  (1) Rhizopus F8E635D0_000401-T1 (annotated Pol X, 507 aa, ~41 bp intron): which intron-open penalty (-J) recovers the annotated
#      1,521 nt? Result: -J 31 and 20 give 1,611 nt (the intron is read as coding), -J 10 and 5 give exactly 1,521 nt.
#  (2) Does -J crash? On Dispira (F1C9E2AA, 3,066 contigs) -J 10 aborts with 'mp_extra_cal: Assertion al == r->qe - r->qs failed' (an
#      internal bug of this miniprot version); this loop records which values complete on that genome with the ligase query set.
# Output: $BENCH/Jtest/. Run: sbatch scripts/04t_test_intron_penalty.sh
set -uo pipefail
: "${SCRATCH:?}"
BENCH=$(readlink -f "$HOME/projects/BFD_dashboard/bench/kx2_ku_family")
CLEAN=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes
MP=/opt/linux/rocky/8.x/x86_64/pkgs/miniprot/0.2/bin/miniprot
w="$SCRATCH/jtest"; mkdir -p "$w" "$BENCH/Jtest"
# Test 1: the Rhizopus Pol X locus at J = 31 20 16 13 10 5: total CDS nucleotides of the model (annotated gene = 1,521 nt)
zcat "$CLEAN/GCA_051997145.1_ASM5199714v1.fa.gz" > "$w/rhizopus.fa"
awk '/^>/{p=($1==">ANN|F8E635D0_000401-T1|F8E635D0")} p' "$BENCH/cand/polx_annotated.faa" > "$w/q1.faa"
for J in 31 20 16 13 10 5; do
  "$MP" -t 2 -G 5k -J "$J" --gff "$w/rhizopus.fa" "$w/q1.faa" 2>/dev/null | awk -F'\t' -v J="$J" '$3=="mRNA"{n++} $3=="CDS"{t+=$5-$4+1} END{print "Rhizopus Pol X, J=" J ": models=" n ", CDS nt=" t}'
done
# Test 2: does -J crash on the Dispira genome (3,066 contigs) with the ligase query set?
zcat "$CLEAN/GCA_027577015.1_UCR_DipaRSA1196_1.0.fa.gz" > "$w/dispira.fa"
cat "$BENCH/refs/lig.faa" "$BENCH/cand/lig_annotated.faa" > "$w/q_lig.faa"
for J in 31 25 20 16 13 10; do
  if "$MP" -t 2 -G 5k -J "$J" --gff "$w/dispira.fa" "$w/q_lig.faa" > "$BENCH/Jtest/dispira_lig_J$J.gff" 2> "$w/err.txt"; then
    echo "Dispira ligase queries, J=$J: ok, $(grep -c -P '\tmRNA\t' "$BENCH/Jtest/dispira_lig_J$J.gff") models"
  else
    echo "Dispira ligase queries, J=$J: CRASH ($(grep -c 'Assertion' "$w/err.txt") assertion)"
  fi
done
