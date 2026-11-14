#!/usr/bin/bash -l
#SBATCH -p short -c 4 --mem 8gb --time 2:00:00 -J kx2_ku_tblastn
#SBATCH --array=1-31%8
#SBATCH -o kx2_ku_tblastn.%A_%a.log
# tblastn of the Ku-pathway query set against the assemblies of the worklist genomes (one array task per genome).
# Submit from the BENCH folder:  BENCH=... REPO=... sbatch $REPO/docs/analyses/kx2_ku_pathway/tblastn/run_tblastn.sh
# Restartable: a genome whose hits file exists is skipped. Temporary BLAST databases go to $SCRATCH (node-local, cleared when
# the task ends); only the small hit tables are kept in $BENCH/hits/.
# Search settings are deliberately permissive (e-value 1e-3) so that weak, diverged hits are recorded; the strong/weak call is
# made later in summarize_tblastn.py, not here.
set -euo pipefail
: "${SCRATCH:?SCRATCH is not set; this job must run under SLURM on UCR HPCC}"
: "${REPO:?set REPO to a checkout of the repo}"; : "${BENCH:?set BENCH to the output folder}"
BENCH=$(readlink -f "$BENCH"); T="$REPO/docs/analyses/kx2_ku_pathway/tblastn"
# Assemblies: input/dna/*.scaffolds.fa are symlinks into genome_annotation/<name>/predict_results/, which no longer exist (checked
# 2026-11-11: 400 of 400 sampled links broken). The cleaned genomes are in input_clean_genomes/<ASMID>.fa.gz (unmasked) and
# <ASMID>.masked.fasta.gz; the layout is described in Fungi_BFD_runs/README_DATA_LAYOUT.md (marked UNVERIFIED there). Read-only.
CLEAN=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes
module load ncbi-blast/2.17.0+
# ROUND=1 (default): queries.faa -> hits/ ; ROUND=2: queries_round2.faa -> hits_round2/ (closer-relative BFD proteins).
# Leave-one-out: a query protein that comes from the genome being searched (accession starts with its locustag) is dropped,
# so a positive control cannot be recovered trivially by its own protein.
# ROUND=3: queries_round3.faa -> hits_round3/, genomes from genomes_round3.tsv (128 genomes: submit with --array=1-128%8, the
# #SBATCH line above is for rounds 1 and 2). Column 6 of genomes_round3.tsv is the expected scaffold count (bfd.asm_stats.contig_count);
# the task fails if the assembly does not match it.
case "${ROUND:-1}" in 1) SUF="";; 2) SUF="_round2";; 3) SUF="_round3";; *) echo "ROUND must be 1, 2 or 3" >&2; exit 2;; esac
GEN="$T/genomes.tsv"; [ "${ROUND:-1}" = 3 ] && GEN="$T/genomes_round3.tsv"
mkdir -p "$BENCH/hits${SUF}"
row=$(tail -n +2 "$GEN" | sed -n "${SLURM_ARRAY_TASK_ID}p")
test -n "$row" || { echo "no genome for task ${SLURM_ARRAY_TASK_ID}" >&2; exit 2; }
loc=$(echo "$row" | cut -f1); name=$(echo "$row" | cut -f2); asmid=$(echo "$row" | cut -f5)
out="$BENCH/hits${SUF}/$loc.tsv"
if [ -s "$out" ]; then echo "$loc already done"; exit 0; fi
asm="$CLEAN/$asmid.fa.gz"
test -s "$asm" || asm="$CLEAN/$asmid.masked.fasta.gz"
test -s "$asm" || { echo "missing assembly for $loc ($asmid) in $CLEAN" >&2; exit 4; }
db="$SCRATCH/db_$loc"; mkdir -p "$db"
zcat "$asm" > "$db/asm.fa"
if [ "${ROUND:-1}" = 3 ]; then
  exp=$(echo "$row" | cut -f6); got=$(grep -c '^>' "$db/asm.fa")
  if [ -n "$exp" ] && [ "$exp" != "$got" ]; then echo "SCAFFOLD COUNT MISMATCH for $loc: $got in file, $exp in asm_stats" >&2; exit 6; fi
fi
makeblastdb -in "$db/asm.fa" -dbtype nucl -out "$db/asm" -parse_seqids >/dev/null
# columns: query, subject scaffold, %identity, aligned length (aa), e-value, bitscore, qstart, qend, qlen, sstart, send, frame
awk -v l="$loc" '/^>/{keep=(index($0, "|" l "_")==0)} keep' "$BENCH/queries${SUF}.faa" > "$db/q.faa"
echo "queries used: $(grep -c '^>' "$db/q.faa") of $(grep -c '^>' "$BENCH/queries${SUF}.faa") (self-genome queries dropped)"
tblastn -query "$db/q.faa" -db "$db/asm" -evalue 1e-3 -max_target_seqs 50 -num_threads 4 \
        -outfmt '6 qseqid sseqid pident length evalue bitscore qstart qend qlen sstart send frames' > "$out.tmp"
mv "$out.tmp" "$out"
echo "$loc ($name, $(basename "$asm")): $(wc -l < "$out") HSPs; assembly $(grep -c '^>' "$db/asm.fa") sequences, $(du -h "$db/asm.fa" | cut -f1)"
