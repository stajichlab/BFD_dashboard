#!/usr/bin/bash -l
#SBATCH -p short -c 4 --mem 16gb --time 2:00:00 -J kickx_final
#SBATCH -o kickx_final.%j.log
# Finalize the improved Kickxellomycota trees (see docs/analyses/kx2_ku_pathway/species_tree/README.md, "What is still running").
# For each dataset: take the PMSF ML tree (pmsf_ml.treefile), add bootstrap support from the PMSF replicates that exist when this runs
# (only if at least MIN_REPS replicates exist; otherwise the tree is reported without PMSF support), then run kickx_tree_report.py.
# The running bootstrap jobs are not touched: their files are copied to a separate directory first.
#
# Usage (from the HPCC bench folder bench/kickx_tree):
#   REPO=~/projects/BFD_dashboard/repo_or_worktree_copy BENCH=~/projects/BFD_dashboard/bench/kickx_tree sbatch --export=ALL kickx_finalize.sh
# Output: $BENCH/iq_final/<mode>/ (copies, supported tree, report) and $BENCH/iq_final/<mode>/report/ (report.txt, rooted Newicks, figures, clade tables).
set -euo pipefail
: "${REPO:?set REPO to a directory holding workflows/phylo/}"
: "${BENCH:?set BENCH to bench/kickx_tree}"
BENCH=$(readlink -f "$BENCH"); REPO=$(readlink -f "$REPO")
MIN_REPS=${MIN_REPS:-50}
module load IQ-TREE/3.0.1
# mode : iq subfolder : taxa table : first-pass tree
declare -A IQ=( [backbone]=backbone [coemansia]=coe_derep )
declare -A TAXA=( [backbone]="$BENCH/taxa/backbone_taxa.tsv" [coemansia]="$BENCH/taxa/coemansia_taxa.tsv" )   # full 121-taxon table: the first-pass tree has all 121 tips; the ML tree uses the 76 dereplicated ones
declare -A FP=(
  [backbone]="$BENCH/backbone/results/protein/buildtree/fungi_odb12/fasttree/protein-kickx_backbone-taxa_63.fungi_odb12.fasttree.support.treefile"
  [coemansia]="$BENCH/coemansia/results/protein/buildtree/fungi_odb12/fasttree/protein-kickx_coemansia-taxa_121.fungi_odb12.fasttree.support.treefile" )
for mode in backbone coemansia; do
  src="$BENCH/iq/${IQ[$mode]}"; dst="$BENCH/iq_final/$mode"
  mkdir -p "$dst"
  cp "$src/pmsf_ml.treefile" "$src/guide.treefile" "$dst/"
  for f in astral.tre astral_c50.tre astral_t2.tre; do [ -f "$src/$f" ] && cp "$src/$f" "$dst/"; done
  if [ -f "$src/guide.contree" ]; then cp "$src/guide.contree" "$dst/"; fi
  reps=0
  if [ -s "$src/pmsf.boottrees" ]; then
    # keep only complete lines (the running job may be mid-write)
    grep -c ';$' "$src/pmsf.boottrees" > "$dst/n_reps.txt" || true
    grep ';$' "$src/pmsf.boottrees" > "$dst/pmsf.boottrees"
    reps=$(cat "$dst/n_reps.txt")
  fi
  echo "$mode: $reps PMSF bootstrap replicates available"
  if [ "$reps" -ge "$MIN_REPS" ]; then
    (cd "$dst" && iqtree3 -sup pmsf_ml.treefile -t pmsf.boottrees --prefix pmsf_ml_sup -redo -quiet)
  else
    echo "$mode: fewer than $MIN_REPS replicates; PMSF support not computed" | tee "$dst/NO_PMSF_SUPPORT.txt"
  fi
  /usr/bin/python3.12 "$REPO/workflows/phylo/kickx_tree_report.py" "$mode" "$dst" "${TAXA[$mode]}" "${FP[$mode]}" "$dst/report"
done
