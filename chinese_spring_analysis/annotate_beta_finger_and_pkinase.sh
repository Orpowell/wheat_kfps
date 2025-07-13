#!/bin/bash/

# Extract Pkinase hmm profile

#hmmfetch /Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/T2T_chinese_spring_analysis/data/Pkinase_hmm/Pfam-A.hmm Pkinase \
#	> /Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/T2T_chinese_spring_analysis/data/Pkinase_hmm/Pkinase.hmm

# Annotate kinase domains in chinese spring non-redundant kfp loci

gene_models="/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/T2T_chinese_spring_analysis/analysis/locus_representatives/T2T_gene_models_and_cloned_kfps.faa"
dir="/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/T2T_chinese_spring_analysis/analysis/KFP_annotation"

hmmsearch \
	--domtblout $dir/T2T_pkinase_domains.txt \
	--cpu 8 \
	--domE 1e-4 \
	/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/T2T_chinese_spring_analysis/data/Pkinase_hmm/Pkinase.hmm \
	$gene_models

# Annotate beta-finger motifs in chinese spring non-redundant kfp loci

hmmsearch \
        --domtblout $dir/T2T_beta_finger_domains.txt \
		--cpu 8 \
		--domE 1e-4 \
        /Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/beta_finger_HMM/analysis/beta_finger_HMM_refinement/refinement/final_output/beta_finger.cs.hmm \
        $gene_models

