
DATADIR="/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/beta_finger_HMM/data/KFP_alphafold3_predictions/folds_2025_06_10_09_00"
OUTDIR="/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/beta_finger_HMM/analysis/structure_quality_assessment"
PTM=$OUTDIR/prediction_ptms.csv

echo "protein,ptm,plddt" > $PTM


for dir in $DATADIR/*/;
do
	protein=$(basename $dir)
	protein_out=$OUTDIR/$protein
	prediction_data=$dir/fold_${protein}_full_data_0.json
	prediction_ptm=$dir/fold_${protein}_summary_confidences_0.json
	prediction_structure=$dir/fold_${protein}_model_0.cif
	
	mkdir -p protein_out
	
	alphafold-analyser pae \
		--alphafold3 \
		--data $prediction_data \
		--output $protein_out/${protein}_pae.png

	alphafold-analyser plddt \
		--alphafold3 \
		--data $prediction_data \
		--output $protein_out/${protein}_plddt.png > ${protein_out}/${protein}_plddt.txt
	
	alphafold-analyser structure \
		--structure $prediction_structure \
		--binary /Users/powellor/anaconda3/envs/alphafold-analyser/bin/pymol \
		--output $protein_out/${protein}_plddt.pse
	
	ptm=$(grep '"ptm"' $prediction_ptm | cut -c 9-12)
	plddt=$(grep "Average pLDDT:" ${protein_out}/${protein}_plddt.txt | cut -c 37-41 )
	
	echo ${protein},${ptm},${plddt} >> $PTM
	
	
done








