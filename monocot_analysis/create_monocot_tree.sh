#!/usr/bin/bash

monocot_species="/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/monocot_KFP_presence/analysis/monocot_genome_species.tsv" 

monocot_tree="/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/monocot_KFP_presence/analysis/monocot_genome_species.nwk"

cut -f3 $monocot_species | ete3 ncbiquery --tree > $monocot_tree


