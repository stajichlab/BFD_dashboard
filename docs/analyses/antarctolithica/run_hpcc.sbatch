#!/bin/bash
#SBATCH --job-name=aoC
#SBATCH --partition=epyc,intel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --output=logs/aoC_%x_%j.log
# Analysis C (Antarctolithica obscura) on the UCR HPCC, one stage per job. Submit all stages, chained, with run_hpcc.sh.
#   sbatch --job-name=aoC_<stage> run_hpcc.sbatch <stage>     stage = prep | rf | phylo | report
# Environment (exported by run_hpcc.sh or set by hand):
#   AO_BASE   working directory (data/, phylo/, mag/, report/, results/, figures); must exist
#   BFD_DATA  shared BFD root (default /bigdata/stajichlab/shared/projects/BFD); read only
#   AO_REPO   repo checkout (required: SLURM runs a spooled copy of this file); scripts = $AO_REPO/docs/analyses/antarctolithica/scripts
set -euo pipefail
STAGE=${1:?stage: prep | rf | phylo | report}
: "${AO_BASE:?set AO_BASE}"
export BFD_DATA=${BFD_DATA:-/bigdata/stajichlab/shared/projects/BFD}
AN=${AO_REPO:?set AO_REPO}/docs/analyses/antarctolithica
S=$AN/scripts
export NCPU=${SLURM_CPUS_PER_TASK:-4}
export PYTHONNOUSERSITE=1   # keep ~/.local site-packages out of the pixi environment
PY="pixi run --manifest-path $AN/pixi.toml python"
cd "$AO_BASE"
echo "stage=$STAGE host=$(hostname) start=$(date -Is) commit=${AO_COMMIT:-unknown} AO_BASE=$AO_BASE BFD_DATA=$BFD_DATA NCPU=$NCPU"
t0=$(date +%s)
step() { local s=$(date +%s); echo "== $* ($(date +%T))"; "$@"; echo "   done in $(( $(date +%s) - s )) s"; }

case $STAGE in
prep)   # data staging, comparison set, features, ordination/PERMANOVA, univariate, SwissProt + contamination, AAI, gene content, composition
    module load diamond/2.1.12
    step $PY $S/00_prep_data.py
    step $PY $S/01_build_set.py                       | tee results/build_set.log
    step $PY $S/01b_swissprot_extract.py
    # phylogeny taxa: re-derive and report, but use the published 159-taxon set (AO_PHYLO_TAXA) so the tree is comparable
    export AO_PHYLO_TAXA=${AO_PHYLO_TAXA:-$AN/results/taxa.tsv}
    step $PY $S/02_phylo_taxa.py                      | tee results/phylo_taxa.log
    cp phylo/taxa.tsv results/taxa.tsv
    step $PY $S/02b_taxa_proteomes.py results/taxa.tsv results/taxa_proteomes_hpcc.tsv
    step $PY $S/04_build_features.py                  | tee results/build_features.log
    step $PY $S/05_ordination.py                      | tee results/ordination.log
    step $PY $S/07_univariate.py                      > results/univariate.log
    step $PY $S/08_swissprot_contam.py                | tee results/swissprot_contam.log
    step $PY $S/08b_mag_aai.py
    step $PY $S/09_gene_content.py                    > results/gene_content.log
    step $PY $S/12_composition.py                     > results/composition.log
    ;;
rf)     # random forests: CV per block, thinning, permutation importance, leave-one-class-out
    step $PY $S/06_rf.py                              > results/rf.log
    ;;
phylo)  # BUSCO odb10 markers -> MAFFT -> trimAl -> concatenation -> FastTree LG+G; jackknife/WAG/reduced-taxon robustness trees
    module load hmmer/3.3.2 mafft/7.505 trimal/1.4.1 fasttree/2.1.11
    export FASTTREE=$(command -v FastTreeMP || command -v FastTree)
    # decompressed proteomes are temporary: keep them on node-local $SCRATCH
    FA=${SCRATCH:-$AO_BASE/phylo}/fa; mkdir -p "$FA"; rm -rf phylo/fa; ln -s "$FA" phylo/fa
    step bash $S/03_phylo_search.sh
    step $PY $S/11_phylo_build.py                     | tee results/phylo_build.log
    cp phylo/concat.fasttree.nwk phylo/marker_presence.csv results/
    step $PY $S/13_phylo_fig.py                       > results/phylo_placement.log
    step $PY $S/14_phylo_jackknife.py
    mkdir -p results/robustness_trees; cp phylo/jk/*.nwk results/robustness_trees/
    step $PY $S/14b_phylo_jk_eval.py
    rm -f phylo/fa
    ;;
report) # figures, report tables, assembled markdown, comparison with the committed laptop results
    step $PY $S/10_figures.py                         > results/figures.log
    step $PY $S/15_rf_fig.py                          > results/rf_fig.log
    step $PY $S/16_report_tables.py
    step $PY $S/17_assemble_report.py
    step $PY $S/18_compare_runs.py $AN/results        > results/rerun_vs_reference.log
    ;;
*) echo "unknown stage $STAGE"; exit 2;;
esac
echo "stage=$STAGE finished $(date -Is) elapsed=$(( $(date +%s) - t0 )) s"
