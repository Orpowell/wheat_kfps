#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1      
#SBATCH -c 16       
#SBATCH -p jic-gpu
#SBATCH --mem 80G
#SBATCH -t 2-12:00
#SBATCH --mail-type=ALL
#SBATCH -o logs/%x/%j/%x.%j.out
#SBATCH -e logs/%x/%j/%x.%j.err
#SBATCH --gres=gpu:1
#SBATCH --array=1-5

# Fixed

land_plants_model="land_plant_v0.3_a_0080.h5"

# Inputs

config="genome_annotation_config.txt"

IFS=',' read -r -a array <<< $(awk -v awkvar="${SLURM_ARRAY_TASK_ID}" 'NR==awkvar' $config)

sample="${array[0]}"
genome="${array[1]}"
out="${array[2]}"
OUTDIR="${out}/${sample}"

mkdir -p $OUTDIR

echo '###INPUTS###'
echo $sample" "$genome" "$out 

gff=$OUTDIR/$sample.gff
proteome=$OUTDIR/$sample.faa
transcriptome=$OUTDIR/$sample.fna

python3 Helixer.py \
	--fasta-path $genome \
	--model-filepath $land_plants_model \
	--subsequence-length 64152 \
	--gff-output-path $gff \
	--batch-size 80

gffread \
	-x $transcriptome \
	-y $proteome \
	-g $genome \
	$gff
