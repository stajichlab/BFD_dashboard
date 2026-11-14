#!/usr/bin/bash -l
#SBATCH -p short -c 4 --mem 8gb --time 2:00:00 -J kx2mt_tblastn
#SBATCH --array=1-128%10
#SBATCH -o kx2mt_tblastn.%A_%a.log
# tblastn of the meiotic-toolkit query set against the assemblies of the 128 genomes (one array task per genome); adapted from
# ../../tblastn/run_tblastn.sh. Submit from the BENCH folder:
#   BENCH=... T_DIR=<folder with genomes.tsv> sbatch --export=ALL run_tblastn.sh
# Restartable (a genome whose hits file exists is skipped). Assemblies are input_clean_genomes/<ASMID>.fa.gz (input/dna links are broken);
# the task fails if the scaffold count differs from bfd.asm_stats.contig_count (column 6 of genomes.tsv). Databases go to $SCRATCH.
# Leave-one-out: a query protein from the genome being searched (accession starts with its locustag) is dropped.
# Search is permissive (e-value 1e-3); strong / weak / none calls are made in summarize_tblastn.py, not here.
set -euo pipefail
: "${SCRATCH:?SCRATCH is not set; this job must run under SLURM on UCR HPCC}"
: "${BENCH:?set BENCH}"; : "${T_DIR:?set T_DIR to the folder that holds genomes.tsv}"
BENCH=$(readlink -f "$BENCH")
CLEAN=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes   # read-only
module load ncbi-blast/2.17.0+
mkdir -p "$BENCH/hits"
row=$(tail -n +2 "$T_DIR/genomes.tsv" | sed -n "${SLURM_ARRAY_TASK_ID}p")
test -n "$row" || { echo "no genome for task ${SLURM_ARRAY_TASK_ID}" >&2; exit 2; }
loc=$(echo "$row" | cut -f1); name=$(echo "$row" | cut -f2); asmid=$(echo "$row" | cut -f5); exp=$(echo "$row" | cut -f6)
out="$BENCH/hits/$loc.tsv"
if [ -s "$out" ]; then echo "$loc already done"; exit 0; fi
asm="$CLEAN/$asmid.fa.gz"; test -s "$asm" || asm="$CLEAN/$asmid.masked.fasta.gz"
test -s "$asm" || { echo "missing assembly for $loc ($asmid)" >&2; exit 4; }
db="$SCRATCH/db_$loc"; mkdir -p "$db"
zcat "$asm" > "$db/asm.fa"
got=$(grep -c '^>' "$db/asm.fa")
if [ -n "$exp" ] && [ "$exp" != "NA" ] && [ "$exp" != "$got" ]; then echo "SCAFFOLD COUNT MISMATCH for $loc: $got in file, $exp in asm_stats" >&2; exit 6; fi
makeblastdb -in "$db/asm.fa" -dbtype nucl -out "$db/asm" -parse_seqids >/dev/null
awk -v l="$loc" '/^>/{keep=(index($0, "|" l "_")==0)} keep' "$BENCH/queries.faa" > "$db/q.faa"
echo "queries used: $(grep -c '^>' "$db/q.faa") of $(grep -c '^>' "$BENCH/queries.faa") (self-genome queries dropped)"
# columns: query, subject scaffold, %identity, aligned length (aa), e-value, bitscore, qstart, qend, qlen, sstart, send, frames
tblastn -query "$db/q.faa" -db "$db/asm" -evalue 1e-3 -max_target_seqs 50 -num_threads 4 \
        -outfmt '6 qseqid sseqid pident length evalue bitscore qstart qend qlen sstart send frames' > "$out.tmp"
mv "$out.tmp" "$out"
echo "$loc ($name, $(basename "$asm")): $(wc -l < "$out") HSPs; $got scaffolds"
