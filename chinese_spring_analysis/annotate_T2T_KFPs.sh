input="/ibex/project/c2210/powellor/kfp_paper/kfp_annotation/data/T2T_KFP/T2T_KFP_full_length.faa"
interpro_data='/ibex/project/c2210/powellor/kfp_paper/kfp_annotation/data/interproscan-5.74-105.0/data'
interpro_sif='/ibex/project/c2210/powellor/kfp_paper/kfp_annotation/data/interproscan_editable.sif'
analysis_dir="/ibex/project/c2210/powellor/kfp_paper/kfp_annotation/analysis/kfp_annotation"

singularity exec \
    -B $interpro_data:/opt/interproscan/data \
    -B $input:/input/T2T_KFP_full_length.faa \
    -B $analysis_dir/temp:/temp \
    -B $analysis_dir/output:/output \
    $interpro_sif \
    /opt/interproscan/interproscan.sh \
    --input /input/T2T_KFP_full_length.faa \
    --disable-precalc \
    --output-dir /output \
    --tempdir /temp \
    --cpu 80
