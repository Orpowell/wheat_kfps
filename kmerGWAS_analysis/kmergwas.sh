#!/bin/bash -e

kmer_pa_matrix_pa_pheno_score_faster \
	-k /jic/common/kmer_matrices/watkins/unified_list.txt \
	-s 31 \
	-p Br48_pwt3_no_rwt3_median_phenotypes.txt \
	-i matrix_list.txt \
	-b 1000 \
	-t 64 \
	-N 1000000 \
	-pop watkins_PCA_full.tsv \
	-npc 9 \
	> kmergwas.txt 
