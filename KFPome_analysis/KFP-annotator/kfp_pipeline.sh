#!/usr/bin/bash

# Variables
SPECIES=$1 # "rice_agis1"
PROTEINS=$2 # "rice_agis1.faa"
OUTDIR=$3 # "OUTDIR"

# Make outdir
mkdir -p $OUTDIR/$SPECIES

BASENAME=$OUTDIR/$SPECIES/$SPECIES

mkdir -p $BASENAME/temp

# HMM models
HMM_BF="beta_finger.monocot.hmm"
HMM_PKINASE="Pkinase.hmm"

#tools
module load hmmer

extract_kfps="python3 extract_KFPs.py"

count_beta_fingers="python3 count_beta_fingers.py"

annotate_kfps="python3 annotate_auxiliary_modules.py"

itol_aux_domains="python3 generate_aux_module_itol.py"

#interpro
interpro_data="interproscan-5.76-107.0/data"
interpro_sif="interproscan.sif"

echo Running HMMsearch...

# Search kinase domains
hmmsearch \
	--domtblout ${BASENAME}_pkinase_domains.txt \
	--cpu 8 \
	--domE 1e-4 \
	$HMM_PKINASE \
	$PROTEINS

# Search beta-finger motifs
hmmsearch \
	--domtblout ${BASENAME}_beta_finger_domains.txt \
	--cpu 8 \
	--domE 1e-4 \
	$HMM_BF \
	$PROTEINS

echo Extracting KFPs...

# Extract KFPs
$extract_kfps \
    ${BASENAME}_beta_finger_domains.txt \
    ${BASENAME}_pkinase_domains.txt \
    $PROTEINS \
    $BASENAME

echo Counting beta-fingers

$count_beta_fingers \
	${BASENAME}_beta_finger_domains.txt \
	${BASENAME}_KFP_full_length.faa \
	$BASENAME

# Build tree

echo Running MUSCLE5...

muscle5 \
	-align ${BASENAME}_KFP_primary_modules.faa \
	-output ${BASENAME}_KFP_primary_modules.aln

#echo Running FastTree...

fasttree \
	${BASENAME}_KFP_primary_modules.aln \
	> ${BASENAME}_KFP_primary_modules.nwk

# Annotate KFPs

echo Running InterProScan...

interproscan.sh \
	--input ${BASENAME}_KFP_full_length.faa \
	--outfile ${BASENAME}_KFP_full_length_annotation.tsv \
	--tempdir $BASENAME/temp \
    --formats TSV \
	--disable-precalc \
	--cpu 64

rm -fr $BASENAME/temp

$annotate_kfps ${BASENAME}_KFP_full_length_annotation.tsv ${BASENAME}_primary_modules.tsv $BASENAME

$itol_aux_domains ${BASENAME}_KFP_auxiliary_modules.tsv $BASENAME