# HMM Refinement

## Contents

### donwload_kfp_sequences.py

Download KFP sequences from NCBI using genbank codes provided.

### generate_alphafold3_prediction_json.py

Generate an AlphaFold3 submission json for each KFP sequence. 

### extract_KFP_beta_fingers.py 

Extract beta-finger sequences based on coordinates provided.

### beta_finger_hmm_auto_refiner.sh 

Iteratively refine a seed HMM on the proteome of Chinese Spring.

### plot_refinement.py  

Plot iterative refinement of seed HMM. (Fig. S1.)

### plot_sequence_logo.py 

Plot sequence logo of HMM sequence from skylign matrix. (Fig. 1A)
