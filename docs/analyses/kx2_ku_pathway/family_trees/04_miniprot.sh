#!/usr/bin/bash -l
#SBATCH -p short -c 4 --mem 12gb --time 2:00:00 -J kx2ft_mp
#SBATCH --array=1-31%8
#SBATCH -o kx2ft_mp.%A_%a.log
# Step 4: rescue gene models for loci that annotation missed. One array task per genome (the 31 worklist genomes in
# strain_genomes.tsv, a copy of ../tblastn/genomes.tsv): map the reference proteins of each family to the assembly with miniprot, keep
# the best model per locus, translate it (04b_translate_models.py). Submit from anywhere:  sbatch scripts/04_miniprot.sh
# Assemblies are input_clean_genomes/<ASMID>.fa.gz (the input/dna links are broken); read-only. Temporary files go to $SCRATCH.
# miniprot 0.2-r116. Maximum intron size is set to 5 kb (-G 5k). The default (200 kb) produced chimeric models that stitched distant
# exons together (first run: ligase models of 92-101 kb in Coemansia, where introns are typically under 400 bp, and 20 of 30 ligase
# loci with internal stops); that run is kept in miniprot_v1_G200k/ and not used.
set -euo pipefail
: "${SCRATCH:?SCRATCH is not set; this job must run under SLURM on UCR HPCC}"
BENCH=$(readlink -f "$HOME/projects/BFD_dashboard/bench/kx2_ku_family")
CLEAN=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes
module load miniprot/0.2 || true
MP=$(command -v miniprot || echo /opt/linux/rocky/8.x/x86_64/pkgs/miniprot/0.2/bin/miniprot)
row=$(tail -n +2 "$BENCH/scripts/strain_genomes.tsv" | sed -n "${SLURM_ARRAY_TASK_ID}p")
test -n "$row" || { echo "no genome for task ${SLURM_ARRAY_TASK_ID}" >&2; exit 2; }
loc=$(echo "$row" | cut -f1); asmid=$(echo "$row" | cut -f5)
mkdir -p "$BENCH/miniprot"
if [ -e "$BENCH/miniprot/$loc.done" ]; then echo "$loc already done"; exit 0; fi
asm="$CLEAN/$asmid.fa.gz"; test -s "$asm" || { echo "missing assembly for $loc ($asmid)" >&2; exit 4; }
w="$SCRATCH/mp_$loc"; mkdir -p "$w"
zcat "$asm" > "$w/asm.fa"
nfail=0
for fam in lig reca polx; do
  # query set = the UniProt references plus every annotated candidate of the family from the BFD genomes (step 3). The distant yeast and
  # animal references alone found no Pol X model at all, even in Rhizopus, where annotation has Pol X proteins; closer fungal relatives
  # make diverged copies findable. A model that re-finds an annotated gene is removed as a duplicate in step 5.
  cat "$BENCH/refs/$fam.faa" "$BENCH/cand/${fam}_annotated.faa" > "$w/q_$fam.faa"
  # -J 16: intron open penalty lowered from 31. These genomes (Mucorales, Kickxellales) have many short introns (about 40-100 bp); with the
  # default penalty miniprot reads them as coding sequence, giving extra residues or internal stops. On the annotated Rhizopus Pol X gene
  # F8E635D0_000401-T1 (507 aa, a ~41 bp intron) -J 31 and 20 give a 1,611 nt model, -J 16, 13, 10 and 5 give exactly the annotated 1,521 nt.
  # -J 10 crashes miniprot 0.2-r116 (assertion in mp_extra_cal) on 13 of 31 genomes, -J 13 and above ran on the crashing test genome; 16 is
  # used for margin (04t_test_intron_penalty.sh).
  # A crash in one family is recorded (<loc>.<fam>.failed) and the other families still run; the genome then gets no .done marker and is
  # listed as 'rescue failed' by step 5 instead of being treated as having no model.
  if "$MP" -t 4 -G 5k -J 16 --gff "$w/asm.fa" "$w/q_$fam.faa" > "$BENCH/miniprot/$loc.$fam.gff" 2> "$w/mp_$fam.err"; then
    /usr/bin/python3.12 "$BENCH/scripts/04b_translate_models.py" "$w/asm.fa" "$loc" "$fam" "$BENCH/miniprot/$loc.$fam.gff" "$BENCH/miniprot/$loc.$fam"
  else
    tail -3 "$w/mp_$fam.err" > "$BENCH/miniprot/$loc.$fam.failed"; rm -f "$BENCH/miniprot/$loc.$fam.gff"; nfail=$((nfail+1))
  fi
done
if [ "$nfail" -eq 0 ]; then touch "$BENCH/miniprot/$loc.done"; else echo "$loc: $nfail family run(s) failed" >&2; fi
