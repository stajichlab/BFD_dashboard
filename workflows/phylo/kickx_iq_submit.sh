#!/usr/bin/bash
# Submit the whole improved-tree pipeline for one taxa set with dependencies (run on the HPCC login node; submits only).
# Usage: NAME=coe_derep RUN=<kickx_tree>/coemansia/results/protein TAXA=<taxa file> REPO=<repo copy> kickx_iq_submit.sh
# Resource requests can be overridden with PMSF_MEM (default 400gb) and PMSF_TIME (default 48:00:00).
set -euo pipefail
: "${NAME:?}"; : "${RUN:?}"; : "${TAXA:?}"; : "${REPO:?}"
BENCH=$(readlink -f "${BENCH:-$HOME/projects/BFD_dashboard/bench/kickx_tree}")
S=$REPO/workflows/phylo/kickx_iqtree.sh
cd "$BENCH/iq"; mkdir -p "$NAME"; cd "$NAME"
EXP=ALL,NAME=$NAME,RUN=$RUN,TAXA=$TAXA,REPO=$REPO,BENCH=$BENCH
sub=$(sbatch --parsable -p short -c 2 --mem 8gb --time 0:30:00 --export=$EXP -J "iq_subset_$NAME" -o "iq_subset.%j.log" "$S" subset)
n_markers_guess=82
guide=$(sbatch --parsable --dependency=afterok:$sub -c 32 --mem 64gb --time 12:00:00 --export=$EXP -J "iq_guide_$NAME" -o "iq_guide.%j.log" "$S" guide)
pmsf=$(sbatch --parsable --dependency=afterok:$guide -c 64 --mem "${PMSF_MEM:-400gb}" --time "${PMSF_TIME:-48:00:00}" --export=$EXP -J "iq_pmsf_$NAME" -o "iq_pmsf.%j.log" "$S" pmsf)
gene=$(sbatch --parsable --dependency=afterok:$sub -p short -c 2 --mem 8gb --time 2:00:00 --array=1-$n_markers_guess --export=$EXP -J "iq_gene_$NAME" -o "iq_gene.%A_%a.log" "$S" gene)
astral=$(sbatch --parsable --dependency=afterany:$gene -p short -c 8 --mem 16gb --time 1:30:00 --export=$EXP -J "iq_astral_$NAME" -o "iq_astral.%j.log" "$S" astral)
echo "$NAME: subset $sub guide $guide pmsf $pmsf gene $gene astral $astral"
