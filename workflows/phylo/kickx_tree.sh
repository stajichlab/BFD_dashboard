#!/usr/bin/bash -l
#SBATCH -p stajichlab -c 2 --mem 8gb --time 24:00:00 -J kickx_tree
#SBATCH -o kickx_tree.%j.log
# Kickxellomycota phylogeny, one tree per call:  sbatch kickx_tree.sh backbone   |   sbatch kickx_tree.sh coemansia
# Steps: (1) select taxa with bin/phylo/kickx_tree_taxa.sql, (2) decompress the chosen proteomes to plain FASTA
# (phyling rejects gzip), (3) run nf_phyling (HMM markers, alignment, concatenation, FastTree).
# This is the FIRST-PASS tree. The site-heterogeneous IQ-TREE and gene-tree summary runs that Analysis E calls for
# take the phyling alignments from this run; see workflows/phylo/README.md.
#
# Everything is written under $BENCH (a purgeable one-off folder). Nothing is written to the shared group folder.
# The nextflow work dir is in $BENCH, not $SCRATCH: $SCRATCH is node-local and cleared when the driver job ends, so it
# cannot hold the shared work dir that -resume needs. Delete $BENCH/<tree>/work when the tree is final.
# Restartable: re-submit the same command; nextflow uses -resume and existing outputs are skipped.
#
# Required:  REPO   a checkout of this repo on the HPCC that contains bin/atlas/run_duckdb_job.py and bin/phylo/
# Optional:  BENCH  default ~/projects/BFD_dashboard/bench/kickx_tree
#            MARKERSET default fungi_odb12 (chosen 2026-11-11; matches BFD's BUSCO lineage). nf_phyling accepts a
#                      comma-separated list, e.g. fungi_odb12,mucoromycota_odb12, but the root taxa (chytrids etc.)
#                      are outside Mucoromycota, so the first pass uses fungi_odb12 only.
# The marker database is downloaded by the PHYLING_DOWNLOAD step into $BENCH/phyling_db, not into the shared atlas
# phyling_db (which holds fungi_odb10 only and is a shared folder we do not write to).
set -euo pipefail
TREE=${1:?usage: sbatch kickx_tree.sh taxa|backbone|coemansia}
case "$TREE" in taxa|backbone|coemansia) ;; *) echo "arg must be taxa, backbone or coemansia" >&2; exit 2;; esac
# Run order: first 'taxa' (selection SQL writes BOTH taxa files, so it must run once, never concurrently), then the
# two tree jobs, which depend on it:
#   T=$(sbatch --parsable -p short --time=0:30:00 kickx_tree.sh taxa)
#   sbatch --dependency=afterok:$T kickx_tree.sh backbone; sbatch --dependency=afterok:$T kickx_tree.sh coemansia
: "${SCRATCH:?SCRATCH is not set; this job must run under SLURM on UCR HPCC}"
: "${REPO:?set REPO to a checkout of the BFD_dashboard repo}"
BENCH=${BENCH:-$HOME/projects/BFD_dashboard/bench/kickx_tree}
MARKERSET=${MARKERSET:-fungi_odb12}
TABLES=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables
# Read-only input from the existing atlas phylo staging: locustag -> proteome path.
PHYSRC=/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/atlas_data/phylo
export NXF_SINGULARITY_CACHEDIR=/bigdata/stajichlab/shared/singularity_cache
export NXF_HOME=$HOME/.nextflow
mkdir -p "$BENCH"
# ~/projects is a symlink into /bigdata on UCR HPCC and the containers only see the real /bigdata path; a
# /rhome/... path appears read-only inside Singularity (PHYLING_DOWNLOAD failed with mkdir: Read-only file system).
BENCH=$(readlink -f "$BENCH")

# (1) taxa. run_duckdb_job.py writes <atlas>/taxa/*.tsv and skips blocks whose output exists.
if [ "$TREE" = taxa ]; then
  /usr/bin/python3.12 "$REPO/bin/atlas/run_duckdb_job.py" \
    --sql "$REPO/bin/phylo/kickx_tree_taxa.sql" --tables "$TABLES" --atlas "$BENCH" \
    --tmp "$SCRATCH/duckdb_tmp" --threads 2 --memory 6GB
  wc -l "$BENCH"/taxa/*.tsv
  exit 0
fi
mkdir -p "$BENCH/$TREE" "$BENCH/phyling_db"
TAXA="$BENCH/taxa/${TREE}_taxa.tsv"
test -s "$TAXA" || { echo "missing $TAXA: run 'kickx_tree.sh taxa' first" >&2; exit 4; }

# (2) proteomes -> plain FASTA named <LOCUSTAG>.fa (the taxon label in the tree is the LOCUSTAG).
PEP="$BENCH/$TREE/pep_fa"; mkdir -p "$PEP"
n=0
while IFS=$'\t' read -r loc path; do
  out="$PEP/$loc.fa"
  if [ ! -s "$out" ]; then zcat "$path" > "$out.tmp" && mv "$out.tmp" "$out"; n=$((n+1)); fi
done < <(awk -F'\t' 'NR==FNR{if(FNR>1) p[$1]=$3; next} FNR>1 && ($2 in p){print $2"\t"p[$2]}' \
           "$PHYSRC/locustag_to_pep.tsv" "$TAXA")
want=$(( $(wc -l < "$TAXA") - 1 )); have=$(ls "$PEP" | wc -l)
echo "staged $n new files; $have of $want taxa have a proteome"
if [ "$have" -ne "$want" ]; then echo "taxa without a proteome path:" >&2
  comm -23 <(cut -f2 "$TAXA" | tail -n +2 | sort) <(ls "$PEP" | sed 's/\.fa$//' | sort) >&2; exit 3; fi

# (3) nf_phyling. Small run (about 65 or 120 proteomes), so far less than the 1,654-genus atlas run that failed.
# Pinned to the same nf_phyling commit as the atlas run. Only flags already used by that run are passed.
cat > "$BENCH/$TREE/phyling.config" <<'EOF'
// Nextflow polls squeue for running jobs. On the busy UCR HPCC that query sometimes returns an incomplete list, so a
// RUNNING job looks like it left the queue; with the default exitReadTimeout (270 s) the task is then failed with
// "terminated for an unknown reason" while the SLURM job keeps running (seen on 2026-11-09 and 2026-11-11). Waiting for
// the .exitcode file instead is safe: it appears when the job really ends.
executor {
    exitReadTimeout = '12 h'
}
process {
    maxRetries = 1
    withName: 'PHYLING_DOWNLOAD' { time = 1.h }
    withName: 'PHYLING_ALIGN'    { time = 12.h; memory = 64.GB; cpus = 16; queue = 'epyc' }
    withName: 'PHYLING_FILTER'   { time = 4.h;  queue = 'epyc' }
    withName: 'PHYLING_TREE'     { time = 8.h;  queue = 'epyc' }
    withName: 'PHYKIT_CONCAT'    { time = 2.h }
    withName: 'FASTTREE'         { time = 8.h;  cpus = 16; queue = 'epyc' }
}
params.fasttree_max_cpus = 16
EOF
module load singularity
module load nextflow/26.04.3
cd "$BENCH/$TREE"
nextflow run stajichlab/nf_phyling -r 52b694f356 \
    -profile singularity_slurm,ucr_hpcc \
    -c "$BENCH/$TREE/phyling.config" \
    -work-dir "$BENCH/$TREE/work" \
    --seq_type protein \
    --input "$PEP" \
    --prefix "kickx_${TREE}" \
    --markerset "$MARKERSET" \
    --phyling_db "$BENCH/phyling_db" \
    --outdir "$BENCH/$TREE/results" \
    -resume
