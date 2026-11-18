#!/usr/bin/bash -l
#SBATCH -p epyc,intel -c 4 --mem 32G --time 1:00:00 -J qs_fig
#SBATCH -o qs_fig.%j.log
# Test, then build the Quarto summary figures and tables from the parquet tables of build_summary.sh.
# Run from the bench folder, with a copy of the repo in $BENCH/code:
#   cd ~/projects/BFD_dashboard/bench/quarto_summary && sbatch code/workflows/summary/make_figures.sh
# Output: $BENCH/out/{figures,tables}/{phylum,subphylum,class,global}; copy these to docs/quarto_summary/ to render.
set -euo pipefail
BENCH=${BENCH:-$HOME/projects/BFD_dashboard/bench/quarto_summary}
PY=${PY:-/usr/bin/python3.12}
cd "$BENCH/code"
$PY -m pytest -q -p no:cacheprovider bin/summary/test_make_figures.py
$PY bin/summary/make_figures.py "$BENCH/data" "$BENCH/out"
