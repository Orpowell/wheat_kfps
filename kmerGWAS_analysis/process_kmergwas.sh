#!/bin/bash -e

# Tools
colourise="python3 colorise_bed.py"
plot_manhattan="python3 plot_manhattan.py"

#Params

OUTDIR="outdir"
BASENAME=$OUTDIR/kmergwas
lengths='lengths.txt' # Path to txt file with length of chromosomes in reference genome
kmers="kmergwas.txt" # KmerGWAS output
blastdb="blastdb" # Path to reference genome blastdb 

awk '{print $1"_"$4"\n>"$1}' $kmers > kmergwas.fna

blastn \
    -query kmergwas.fna \
    -db $blastdb \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen" \
    -out ${BASENAME}.blastn \
    -num_threads $SLURM_CPUS_PER_TASK \
    -word_size 31

awk '{split($1, array, "_"); print $2"\t"$9"\t"$10"\t"array[1]"\t"array[2]}' ${BASENAME}.blastn > ${BASENAME}.bed

$colourise ${BASENAME}.bed > ${BASENAME}.coloured.bed

$plot_manhattan $BASENAME.coloured.bed $lengths $BASENAME.manhattan.png 
