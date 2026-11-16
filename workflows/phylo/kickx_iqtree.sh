#!/usr/bin/bash -l
#SBATCH -p epyc -c 32 --mem 400gb --time 48:00:00 -J kickx_iq
#SBATCH -o kickx_iq.%x.%j.log
# Site-heterogeneous ML tree (IQ-TREE 3.0.1, LG+C20+F+G with PMSF) and gene-tree summary (ASTRAL 5.15.5) on the nf_phyling alignments.
# One stage per call, normally submitted with dependencies by kickx_iq_submit.sh:
#   subset  restrict the per-marker alignments of a kickx_tree run to a taxa list; rebuild concat.faa            (short job)
#   guide   LG+F+G ML tree with 1000 UFBoot + SH-aLRT (guide tree for PMSF, and a fast-model comparison tree)
#   pmsf    LG+C60+F+G under the PMSF approximation on the guide tree, 100 nonparametric bootstraps (writes pmsf.boottrees; ~35 min per
#           replicate for 63 taxa x 51,800 sites, so ~2 days; do NOT use pmsf.treefile as the ML tree, see pmsfml)
#   pmsfml  the same PMSF ML tree search without bootstrap (pmsf_ml.treefile is the ML tree)
#   gene    array: one IQ-TREE per marker (ModelFinder among LG,WAG,JTT,Q.PFAM; 1000 UFBoot)
#   astral  ASTRAL 5.15.5 on the gene trees (all branches) and on trees with branches of UFBoot < 50 collapsed (sensitivity)
# Usage: NAME=coe_derep RUN=<kickx_tree>/coemansia/results/protein TAXA=<taxa file> REPO=<repo copy> sbatch [opts] kickx_iqtree.sh <stage>
# Everything is written to $BENCH/iq/$NAME (purgeable).  Mixture model LG+C60+F+G (override with MODEL=LG+C20+F+G). Memory scales with
# classes: the LG+F+G guide fit reports about 1.6 GB for 76 taxa x 39,700 sites, so C60 needs on the order of 100 GB (160 GB requested);
# an earlier draft of this script assumed 250-400 GB, which was an over-estimate and delayed the jobs in the queue.
set -euo pipefail
STAGE=${1:?stage: subset|guide|pmsf|gene|astral}
: "${NAME:?}"; : "${RUN:?}"; : "${TAXA:?}"; : "${REPO:?}"
BENCH=$(readlink -f "${BENCH:-$HOME/projects/BFD_dashboard/bench/kickx_tree}")
D=$BENCH/iq/$NAME; mkdir -p "$D"
T=${SLURM_CPUS_PER_TASK:-8}
case "$STAGE" in
  subset)
    /usr/bin/python3.12 "$REPO/workflows/phylo/kickx_subset_aln.py" "$RUN" "$TAXA" "$D" 8 ;;
  guide)
    module load IQ-TREE/3.0.1
    cd "$D"
    iqtree3 -s concat.faa -m LG+F+G -B 1000 -alrt 1000 -T "$T" --prefix guide --seed 20261112 -redo
    iqtree3 --version | head -1 > versions.txt ;;
  pmsf)
    module load IQ-TREE/3.0.1
    cd "$D"
    iqtree3 -s concat.faa -m "${MODEL:-LG+C60+F+G}" -ft guide.treefile -b 100 -T "$T" --prefix pmsf --seed 20261112 -redo ;;
  pmsfml)
    # ML tree only (no bootstrap). With -b, IQ-TREE 3.0.1 overwrites <prefix>.treefile with the tree of every bootstrap replicate and
    # keeps no copy of the original ML tree, so the ML tree must come from a run without -b. Support for this tree is then mapped from
    # the replicates of the 'pmsf' stage:  iqtree3 -sup pmsf_ml.treefile -t pmsf.boottrees --prefix pmsf_ml_sup
    module load IQ-TREE/3.0.1
    cd "$D"
    iqtree3 -s concat.faa -m "${MODEL:-LG+C60+F+G}" -ft guide.treefile -T "$T" --prefix pmsf_ml --seed 20261112 -redo ;;
  gene)
    module load IQ-TREE/3.0.1
    cd "$D"; mkdir -p gene_trees
    g=$(ls genes | sed -n "${SLURM_ARRAY_TASK_ID}p")
    [ -n "$g" ] || exit 0          # array is sized for the maximum number of markers
    iqtree3 -s "genes/$g" -m MFP -mset LG,WAG,JTT,Q.PFAM -B 1000 -T 2 --prefix "gene_trees/${g%.fa}" --seed 20261112 -redo -quiet ;;
  astral)
    module load ASTRAL/5.15.5
    cd "$D"
    cat gene_trees/*.treefile > genes_all.tre
    /usr/bin/python3.12 "$REPO/workflows/phylo/kickx_tree_tools.py" collapse genes_all.tre genes_c50.tre 50
    # default output: local posterior probability as the branch label; -t 2: full quartet annotation (q1-q3, pp1-pp3, ...)
    java -jar "$ASTRALJAR" -i genes_all.tre -o astral.tre 2> astral.log
    java -jar "$ASTRALJAR" -i genes_all.tre -o astral_t2.tre -t 2 2> astral_t2.log
    java -jar "$ASTRALJAR" -i genes_c50.tre -o astral_c50.tre 2> astral_c50.log
    echo "genes: $(wc -l < genes_all.tre)" >> astral.log ;;
  *) echo "unknown stage $STAGE" >&2; exit 2 ;;
esac
