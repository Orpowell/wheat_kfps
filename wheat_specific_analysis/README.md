# Wheat_specific_analysis

## Contents

### Kinome/extract_wheat_kinome.py

Extract the kinome of chinese spring and generate an iTOL annotation for kinase domains containing an extended beta-finger (Fig. S27)

### Kinome/build_kinome_tree.sh

Build tree of chinese spring kinome (Fig. S27).

### Kinome/blastp_yr15_CS.sh

BLAST kinase domains of YR15 against kinome to identify homologues (Fig. S27). iTOL annotation for kinase domains was generated manually using the best hit. 

### KFP_positions/analyse_kfp_positions.py

Generates a karyogram of the T2T Chinese Spring assembly highlighting the centromeres and
positions of KFPs and iTOL annotation of KFP subgenome association (Fig. S28-29).

### KFP_positions/calculate_nlr_proximity.py

Generates a cumulative density plot of distance of each KFP to nearest CC-NLR gene (Fig. S32).

### pim_analysis/analyse_pim.py

Generate percentage-identity matrices for MSP and Kinase-KFPs (Fig. S30-31). 
