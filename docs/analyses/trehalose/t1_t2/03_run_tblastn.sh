#!/usr/bin/bash -l
#SBATCH -p short -c 4 --mem 8gb --time 2:00:00 -J t1t2_tblastn
#SBATCH --array=1-33%8
#SBATCH -o t1t2_tblastn.%A_%a.log
# tblastn of queries.faa against each assembly in genomes.tsv (one array task per genome). Submit from the bench folder:
#   BENCH=$PWD sbatch --array=1-N%8 code/docs/analyses/trehalose/t1_t2/03_run_tblastn.sh      (N = rows of genomes.tsv)
# Restartable (a genome with a hits file is skipped). Temporary BLAST databases go to $SCRATCH, only the small hit tables are kept.
# Permissive search (e-value 1e-3, up to 50 targets, default masking); the strong/weak/none call is made in 04_summarize.py.
# A query whose accession starts with the locustag of the searched genome is dropped (leave-one-out).
set -euo pipefail
: "${SCRATCH:?SCRATCH is not set; this job must run under SLURM on UCR HPCC}"
: "${BENCH:?set BENCH to the bench folder that holds genomes.tsv and queries.faa}"
BENCH=$(readlink -f "$BENCH")
CLEAN=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes     # read-only
module load ncbi-blast/2.17.0+
mkdir -p "$BENCH/hits"
row=$(tail -n +2 "$BENCH/genomes.tsv" | sed -n "${SLURM_ARRAY_TASK_ID}p")
test -n "$row" || { echo "no genome for task ${SLURM_ARRAY_TASK_ID}" >&2; exit 2; }
loc=$(echo "$row" | cut -f1); asmid=$(echo "$row" | cut -f5)
out="$BENCH/hits/$loc.tsv"
if [ -s "$out" ] || [ -e "$out.done" ]; then echo "$loc already done"; exit 0; fi
asm="$CLEAN/$asmid.fa.gz"; test -s "$asm" || asm="$CLEAN/$asmid.masked.fasta.gz"
test -s "$asm" || { echo "missing assembly for $loc ($asmid)" >&2; exit 4; }
db="$SCRATCH/db_$loc"; mkdir -p "$db"
zcat "$asm" > "$db/asm.fa"
makeblastdb -in "$db/asm.fa" -dbtype nucl -out "$db/asm" -parse_seqids >/dev/null
awk -v l="$loc" '/^>/{keep=(index($0, "|" l "_")==0)} keep' "$BENCH/queries.faa" > "$db/q.faa"
echo "queries used: $(grep -c '^>' "$db/q.faa") of $(grep -c '^>' "$BENCH/queries.faa")"
# columns: query, subject scaffold, %identity, aligned length (aa), e-value, bitscore, qstart, qend, qlen, sstart, send, frames
tblastn -query "$db/q.faa" -db "$db/asm" -evalue 1e-3 -max_target_seqs 50 -num_threads 4 \
        -outfmt '6 qseqid sseqid pident length evalue bitscore qstart qend qlen sstart send frames' > "$out.tmp"
mv "$out.tmp" "$out"; touch "$out.done"
echo "$loc ($asmid): $(wc -l < "$out") HSPs; $(grep -c '^>' "$db/asm.fa") scaffolds"
