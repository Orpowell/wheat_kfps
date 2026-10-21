# HMM Refinement

## Contents

### donwload_kfp_sequences.py

Download KFP sequences from NCBI using genbank codes provided.

### generate_alphafold3_prediction_json.py

Generate an AlphaFold3 submission json for each KFP sequence. 

### extract_KFP_beta_fingers.py 

Extract beta-finger sequences based on coordinates provided (See Table S2).

### build_seed_hmm.sh

Build Seed HMM from the beta-finger sequences of cloned KFPs.

### assess_structure_quality.sh

Assess the quality of AlphaFold3 predictions of KFPs using AlphaFold-Analyser.

### download_monocot_proteomes.sh

Download all available monocot proteomes from NCBI.

### beta_finger_hmm_auto_refiner.sh 

Iteratively refine a seed HMM on the proteome of Chinese Spring.

### plot_refinement.py  

Plot iterative refinement of seed HMM. (Fig. S21).

### plot_sequence_logo.py 

Plot sequence logo of HMM sequence from skylign matrix. (Fig. 1B)


