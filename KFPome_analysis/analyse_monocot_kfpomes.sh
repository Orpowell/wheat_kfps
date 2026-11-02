#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1      
#SBATCH -c 16       
#SBATCH -p jic-short
#SBATCH --mem 32G 
#SBATCH -t 0-10:00
#SBATCH --mail-type=ALL
#SBATCH -o logs/%x/%j/%x.%j.out
#SBATCH -e logs/%x/%j/%x.%j.err
#SBATCH --array=1-5

config="/jic/scratch/groups/Paul-Nicholson/ollie/rmg2/data/monocot_crops/monocot_proteome_config.txt"
pipeline="/jic/scratch/groups/Paul-Nicholson/ollie/rmg2/scripts/monocot_crops/kfp_pipeline.sh"

IFS=',' read -r -a array <<< $(awk -v awkvar="${SLURM_ARRAY_TASK_ID}" 'NR==awkvar' $config)

sample="${array[0]}" 
proteome="${array[1]}"
out="${array[2]}"

echo '###INPUTS###'
echo $sample
echo $proteome
echo $out 

$pipeline $sample $proteome $out

