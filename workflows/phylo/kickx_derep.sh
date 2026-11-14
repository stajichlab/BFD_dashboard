#!/usr/bin/bash -l
#SBATCH -p short -c 16 --mem 16gb --time 2:00:00 -J kickx_derep
#SBATCH -o kickx_derep.%j.log
# Step 1 of the improved Coemansia tree: ANI dereplication of the 107 Coemansia strains of taxa/coemansia_taxa.tsv.
# fastANI 1.34 all-vs-all (default k-mer 16, fragment length 3000) on the cleaned assemblies, then kickx_derep.py clusters at
# >= 99 % ANI (pre-stated; 97 and 95 % reported as sensitivity). Writes into $BENCH/iq/derep (purgeable).
# Usage (from anywhere, on the HPCC):  REPO=<checkout> sbatch kickx_derep.sh
set -euo pipefail
: "${REPO:?set REPO to a checkout of the BFD_dashboard repo}"
BENCH=$(readlink -f "${BENCH:-$HOME/projects/BFD_dashboard/bench/kickx_tree}")
GEN=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes     # read-only
OUT=$BENCH/iq/derep; mkdir -p "$OUT"
TAXA=$BENCH/taxa/coemansia_taxa.tsv
module load fastani/1.34
awk -F'\t' -v g="$GEN" 'NR>1 && $1=="ingroup"{print g"/"$3".fa.gz"}' "$TAXA" > "$OUT/genomes.list"
n=$(wc -l < "$OUT/genomes.list"); missing=$(while read -r f; do test -s "$f" || echo "$f"; done < "$OUT/genomes.list" | wc -l)
echo "$n genomes listed, $missing missing"; test "$missing" -eq 0
fastANI --ql "$OUT/genomes.list" --rl "$OUT/genomes.list" -t "${SLURM_CPUS_PER_TASK:-16}" -o "$OUT/fastani.out" 2> "$OUT/fastani.err"
fastANI --version 2>&1 | head -1 > "$OUT/versions.txt" || true
/usr/bin/python3.12 "$REPO/workflows/phylo/kickx_derep.py" "$TAXA" "$OUT/fastani.out" "$OUT" 99.0
