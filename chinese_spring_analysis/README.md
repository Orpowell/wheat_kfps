# Chinese Spring Analysis

## Contents

### get_locus_representatives.py

Get longest gene model from each locus in Chinese Spring.

### get_pfam.sh

Download Pfam and extract the Pkinase HMM.

### annotate_beta_finger_and_pkinase.sh

Annotate longest gene models using extended beta-fingers and pkinase HMMs.

### extract_kfp_sequences.py

Extract sequences containing an extended beta-finger kinase domain from gene model representatives.

### plot_KFP_physical_positions.py 

Plot physical positions of each KFP in Chinese Spring. (Fig. 1C and Table S4-5)

### annotate_subgenome.py

Generate an iTOL annotation of the subgenome of origin of each KFP. (Fig. S22)

### annotate_T2T_KFPs.sh

Annotate all KFP sequences using InterProScan including TMHMM.

### tree_construction.sh

Generate a phylogenetic tree of Chinese Spring KFPs using only their extended beta-finger kinase domain(s). (Fig. 1D)

### annotate_KFP_auxiliary_domain.py

Extract auxiliary module of each KFP and generate an iTOL annotation. (Fig. 1D, Supplementary Fig. 23A & 24A and Table S6-7)

### annotate_cloned_kfps.py

Generate an iTOL annotation of cloned KFPs. (Fig. 1D)

### count_beta_fingers.py

Count extened beta-finger kinase domains associated with each KFP and generate an itol annotation. (Fig. 1D and Table S8)

### generate_PIM.sh

Generate a percentage identity matrix from a multiple sequence alignment.

### calculate_nlr_proximity.py

Find and plot the proximal CC-NLR to each KFP in chinese spring. (Table S9, Supplementary Fig. 25)

### Analyse_PIM.py

Plot cluster maps of MSP and Kinase KFPs. (Supplementary Fig. 23B & 24B)

