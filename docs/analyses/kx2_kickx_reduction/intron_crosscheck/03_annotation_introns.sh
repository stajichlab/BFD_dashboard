#!/usr/bin/bash -l
#SBATCH -p short -c 4 --mem 32gb --time 1:00:00 -J kx2_ic_annot
#SBATCH -o logs/kx2_ic_annot.%j.log
# Submit from $BENCH:  BENCH=... sbatch $BENCH/scripts/03_annotation_introns.sh
set -euo pipefail
: "${BENCH:?set BENCH}"; BENCH=$(readlink -f "$BENCH"); cd "$BENCH"
/usr/bin/python3.12 scripts/03_annotation_introns.py genomes_resolved.tsv annotation_introns.tsv
