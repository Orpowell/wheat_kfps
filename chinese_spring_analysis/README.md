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

plot physical positions of each KFP in Chinese Spring. (Fig. 1B and Table S2-3)

### annotate_T2T_KFPs.sh

Annotate all KFP sequences using InterProScan including TMHMM.

### tree_construction.sh

Generate a phylogenetic tree of Chinese Spring KFPs using only their extended beta-finger kinase domain(s). (Fig. 1C)

### annotate_KFP_auxiliary_domain.py

Extract auxiliary module of each KFP and generate an iTOL annotation. (Fig. 1C, Supplementary Fig. 3-4 and Table S4-5)

### annotate_cloned_kfps.py

Generate an iTOL annotation of cloned KFPs. (Fig. 1C)

### annotate_subgenome.py

Generate an iTOL annotation of the subgenome of each KFP. (Supplementary Fig. 2)

### count_beta_fingers.py

Count extened beta-finger kinase domains associated with each KFP and generate an itol annotation. (Fig. 1C and Table S6)

### calculate_nlr_proximity.py

Find and plot the proximal CC-NLR to each KFP in chinese spring. (Table S7, Supplementary Fig. 5)


