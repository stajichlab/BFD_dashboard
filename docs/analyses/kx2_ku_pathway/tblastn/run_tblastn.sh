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
MAP=/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/atlas_data/phylo/locustag_to_pep.tsv   # read-only
DNA=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input/dna                                # read-only
module load ncbi-blast/2.17.0+
mkdir -p "$BENCH/hits"
loc=$(tail -n +2 "$T/genomes.tsv" | sed -n "${SLURM_ARRAY_TASK_ID}p" | cut -f1)
test -n "$loc" || { echo "no genome for task ${SLURM_ARRAY_TASK_ID}" >&2; exit 2; }
out="$BENCH/hits/$loc.tsv"
if [ -s "$out" ]; then echo "$loc already done"; exit 0; fi
name=$(awk -F'\t' -v l="$loc" '$1==l{print $2}' "$MAP"); test -n "$name" || { echo "$loc not in locustag map" >&2; exit 3; }
asm="$DNA/$name.scaffolds.fa"; test -s "$asm" || { echo "missing assembly $asm" >&2; exit 4; }
db="$SCRATCH/db_$loc"; mkdir -p "$db"
makeblastdb -in "$asm" -dbtype nucl -out "$db/asm" -parse_seqids >/dev/null
# columns: query, subject scaffold, %identity, aligned length (aa), e-value, bitscore, qstart, qend, qlen, sstart, send, frame
tblastn -query "$BENCH/queries.faa" -db "$db/asm" -evalue 1e-3 -max_target_seqs 50 -num_threads 4 \
        -outfmt '6 qseqid sseqid pident length evalue bitscore qstart qend qlen sstart send frames' > "$out.tmp"
mv "$out.tmp" "$out"
echo "$loc ($name): $(wc -l < "$out") HSPs; assembly $(grep -c '^>' "$asm") sequences, $(du -h "$asm" | cut -f1)"
