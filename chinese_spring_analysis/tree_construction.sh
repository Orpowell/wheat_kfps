muscle \
	-align /ibex/project/c2210/powellor/kfp_paper/kfp_phylogenetics/data/T2T_KFP_primary_modules.faa \
	-output /ibex/project/c2210/powellor/kfp_paper/kfp_phylogenetics/analysis/T2T_KFP_primary_modules.aln \
	-threads 60

fasttree /ibex/project/c2210/powellor/kfp_paper/kfp_phylogenetics/analysis/T2T_KFP_primary_modules.aln > /ibex/project/c2210/powellor/kfp_paper/kfp_phylogenetics/analysis/T2T_KFP_primary_modules.nwk

