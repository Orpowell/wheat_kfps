#!/usr/bin/bash

clustalo \
	-i /Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/T2T_chinese_spring_analysis/analysis/KFP_phylogenetics/T2T_KFP_primary_modules.aln \
	--percent-id \
	--distmat-out=/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/T2T_chinese_spring_analysis/analysis/KFP_phylogenetics/T2T_KFP_primary_modules.pim \
	--full \
	--threads 6 \
	--force
