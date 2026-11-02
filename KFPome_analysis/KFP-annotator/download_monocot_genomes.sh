#!/bin/bash -e
#SBATCH -p nbi-download
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem 4G
#SBATCH -t 0-06:00
#SBATCH --output=logs/%x/%j/%x.%j.out
#SBATCH --error=logs/%x/%j/%x.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ollie.powell2@jic.ac.uk

source activate /hpc-home/yok25vun/mamba/envs/ncbi_datasets

OUTDIR="/jic/scratch/groups/Paul-Nicholson/ollie/rmg2/data/monocot_crops/genomes"

datasets download genome accession GCF_902167145.1 GCA_040256815.2 GCA_951802345.1 GCF_904849725.1 GCF_034140825.1 --include genome --filename $OUTDIR/poaceae_crops.zip

