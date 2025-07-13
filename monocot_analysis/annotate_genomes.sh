pfam="/ibex/project/c2210/powellor/kfp_paper/monocot_presence_absence/data/beta_finger.cs.hmm"

genome_data="/ibex/project/c2210/powellor/kfp_paper/monocot_presence_absence/data/genomes.txt"

genome=$(awk -v task=$SLURM_ARRAY_TASK_ID 'NR==task {print $1}' $genome_data)

base=$(basename $genome .fna)

outdir="/ibex/project/c2210/powellor/kfp_paper/monocot_presence_absence/analysis/annotations/${base}"

mkdir -p $outdir

cd outdir

hmm-annotator \
    --pfam $pfam \
    --sequences $genome \
    --cores 10 \
    --output $outdir/$base.txt \
    --bed $outdir/$base.bed
