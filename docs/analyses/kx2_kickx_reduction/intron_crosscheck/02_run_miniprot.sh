#!/usr/bin/bash -l
#SBATCH -p short -c 4 --mem 16gb --time 2:00:00 -J kx2_ic_miniprot
#SBATCH --array=1-83%10
#SBATCH -o logs/kx2_ic.%A_%a.log
# miniprot (spliced protein-to-genome alignment, no gene finder) of the fixed BUSCO fungi_odb12 ancestral protein set against
# one assembly per array task, with three splice settings, then intron calls with miniprot_introns.py.
# Run from $BENCH (a purgeable folder):  BENCH=... sbatch $BENCH/scripts/02_run_miniprot.sh      (mkdir -p $BENCH/logs first)
# Restartable: a genome whose three tables exist is skipped. Temp files are in node-local $SCRATCH.
# Settings: A default (-J 31 -C 11 -G 200k); B permissive (-J 16 -C 6 -G 5000); C very permissive (-J 8 -C 3 -G 5000).
set -euo pipefail
: "${SCRATCH:?SCRATCH is not set; this job must run under SLURM on UCR HPCC}"
: "${BENCH:?set BENCH}"; BENCH=$(readlink -f "$BENCH"); S="$BENCH/scripts"
module load miniprot/0.2
row=$(tail -n +2 "$BENCH/genomes_resolved.tsv" | sed -n "${SLURM_ARRAY_TASK_ID}p")
test -n "$row" || { echo "no genome for task ${SLURM_ARRAY_TASK_ID}" >&2; exit 2; }
loc=$(echo "$row" | cut -f1); genus=$(echo "$row" | cut -f2); contigs=$(echo "$row" | cut -f8); asm=$(echo "$row" | cut -f12)
mkdir -p "$BENCH/gff" "$BENCH/tables"
if [ -s "$BENCH/tables/A_$loc.tsv" ] && [ -s "$BENCH/tables/B_$loc.tsv" ] && [ -s "$BENCH/tables/C_$loc.tsv" ]; then echo "$loc done"; exit 0; fi
w="$SCRATCH/ic_$loc"; mkdir -p "$w"
zcat "$asm" > "$w/asm.fa"
got=$(grep -c '^>' "$w/asm.fa")
if [ "$got" != "$contigs" ]; then echo "SCAFFOLD COUNT MISMATCH for $loc: $got in file, $contigs in asm_stats" >&2; exit 6; fi
run() {  # setting name, miniprot options
  local name=$1; shift
  local g="$BENCH/gff/${name}_$loc.gff.gz"
  # the GFF is written under a temporary name and renamed only when miniprot finished (pipefail), so a cancelled task can
  # never leave a truncated file that a later run would reuse; an existing complete GFF is reused and only re-parsed
  if [ ! -s "$g" ]; then
    miniprot -t 4 --gff "$@" "$w/asm.fa" "$BENCH/ref/ancestral.faa" 2> "$w/mp_$name.err" | gzip > "$g.tmp"
    mv "$g.tmp" "$g"
  fi
  python3 "$S/miniprot_introns.py" --gff "$g" --ref "$BENCH/ref/ancestral.faa" --genome "$w/asm.fa" \
          --out "$BENCH/tables/${name}_$loc.tsv.tmp" 2>> "$w/mp_$name.err"
  mv "$BENCH/tables/${name}_$loc.tsv.tmp" "$BENCH/tables/${name}_$loc.tsv"
}
run A
run B -J 16 -C 6 -G 5000
# miniprot aborts (core dump) with very low intron penalties on some assemblies (seen for Cryptocalciella; fork A saw the same with
# -J 10), so a failure of setting C is recorded and that genome is left out of C-based comparisons, instead of failing the task
if ! run C -J 8 -C 3 -G 5000; then echo -e "$loc\t$genus\tC_failed" >> "$BENCH/failures.tsv"; rm -f "$BENCH/gff/C_$loc.gff.gz.tmp"; fi
nl() { [ -s "$1" ] && echo $(( $(wc -l < "$1") - 1 )) || echo NA; }
echo "$loc ($genus): scaffolds $got ok; proteins aligned A/B/C: $(nl "$BENCH/tables/A_$loc.tsv") / $(nl "$BENCH/tables/B_$loc.tsv") / $(nl "$BENCH/tables/C_$loc.tsv")"
