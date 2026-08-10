#!/usr/bin/bash -l
#SBATCH -p batch -c 2 --mem 4gb --out logs/noinvenio.%A.log

module load nextflow
module load singularity
export NXF_SINGULARITY_CACHE=/bigdata/stajichlab/shared/singularity_cache/

SRC=/bigdata/stajichlab/jstajich/projects/NovInvenio
nextflow run $SRC/main.nf \
      -profile slurm \
      -c $SRC/conf/ucr_hpcc_slurm.config \
      --config configs/Rmuc.csv \
      --data_dir input \
      --run_tool diamond \
      --cluster_tool pairwise \
      --pfam_hmm db/pfam/Pfam-A.hmm \
      --swissprot_dmnd db/uniprot/uniprot_sprot.fasta.dmnd \
      --modelorgs_config $SRC/configs/modelorgs.yaml \
      --outdir results/NovInvenio -resume

nextflow run $SRC/main.nf \
      -profile slurm \
      -c $SRC/conf/ucr_hpcc_slurm.config \
      --config configs/Rhodotorula.csv \
      --data_dir input \
      --run_tool diamond \
      --cluster_tool pairwise \
      --pfam_hmm db/pfam/Pfam-A.hmm \
      --swissprot_dmnd db/uniprot/uniprot_sprot.fasta.dmnd \
      --modelorgs_config $SRC/configs/modelorgs.yaml \
      --outdir results/NovInvenio -resume
