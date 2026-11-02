#!/bin/bash

blastp \
    -query yr15_kinase_domains.faa \
    -subject wheat_chinese_spring_kinome.faa \
    -num_threads $SLURM_CPUS_PER_TASK \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen" \
    -out Yr15_cs_kinome.blastp
