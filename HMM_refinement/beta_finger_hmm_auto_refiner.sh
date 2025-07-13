#!/bin/bash


current_hmm_profile="/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/beta_finger_HMM/analysis/beta_finger_seed_HMM/seed_beta_finger.hmm"
sequence_db="/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/beta_finger_HMM/data/T2T_chinese_spring/ncbi_dataset/data/GCA_040256815.2/protein.faa"
workdir="/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/beta_finger_HMM/analysis/beta_finger_HMM_refinement/refinement"

rm -r $workdir

mkdir $workdir

prior_hits_1=0
prior_hits_2=0
iteration=1
not_converged=true

while $not_converged
do
	
	# Make iteration work dir
	dir=$workdir/iteration_$iteration
	mkdir -p $dir

	# create file names
	hmmsearch_out=$dir/chinese_spring_refinement.iteration_$iteration.txt
 	
	hmmsearch_aln=$dir/chinese_spring_refinement.iteration_$iteration.aln

 	sequences=$dir/chinese_spring_beta_fingers.iteration_$iteration.faa       
 	
	new_profile=$dir/beta_finger.iteration_$iteration.hmm
	
	hmmsearch --incdomE 1e-10 --cpu 9 -A $hmmsearch_aln --tblout $hmmsearch_out  $current_hmm_profile $sequence_db > $dir/hmmsearch.log
	
	hmmbuild -n beta_finger_iteration_$iteration $new_profile $hmmsearch_aln > $dir/hmmbuild.out 2> $dir/hmmbuild.err
	
	hmmpress $new_profile > $dir/hmmpress.out 2> $dir/hmmpress.err
	
	unique_hits=$(awk 'NR==13{print $3}' $dir/hmmbuild.out)
	
	echo $iteration $unique_hits	

	if [[ $unique_hits -eq $prior_hits_1 ]] &&  [[ $unique_hits -eq $prior_hits_2 ]]; then
		not_converged=false
 	
		mkdir -p $workdir/final_output

		hmmbuild -n beta_finger $workdir/final_output/beta_finger.cs.hmm $hmmsearch_aln > $workdir/final_output/hmmbuild.out 2> $workdir/final_output/hmmbuild.err

        hmmpress $workdir/final_output/beta_finger.cs.hmm > $workdir/final_output/hmmpress.final.out 2> $workdir/final_output/hmmpress.err
	
	else
	 	prior_hits_2=$prior_hits_1
		prior_hits_1=$unique_hits
		iteration=$((iteration + 1))
		current_hmm_profile=$new_profile
	fi
	
done
	






