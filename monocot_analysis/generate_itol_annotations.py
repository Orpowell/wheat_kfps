#!/usr/bin/python3
from glob import glob
import pandas as pd

def count2bool(x):
	if int(x) > 0:
		return 1
	else:
		return 0

species_data = "/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/monocot_KFP_presence/analysis/monocot_genome_species.tsv"
annotation_dir = "/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/monocot_KFP_presence/analysis/annotations"

output_beta_finger_presence="/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/monocot_KFP_presence/analysis/monocot_beta_finger_presence.itol"
output_poaceae_label="/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/monocot_KFP_presence/analysis/poaceae.itol"
output_species_labels="/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/monocot_KFP_presence/analysis/monocot_species_labels.itol"

species_beta_finger_count = dict()

for file in glob(f"{annotation_dir}/*/*.bed"):

	genome = "_".join(file.split("/")[-1].split("_")[:2])
	try:
		data = pd.read_csv(file, sep="\t", header=None)
		data["length"] = data[2] - data[1]
		data = data[(data.length >= 90) & (data[4] <= 1e-10)]
		species_beta_finger_count[genome] = len(data)

	except pd.errors.EmptyDataError:
		species_beta_finger_count[genome] = 0

species_df = pd.read_csv(species_data, sep="\t", header=None)
species_df.columns =["species", "genome", "taxid"]

species_df["bf_count"] = species_df["genome"].map(species_beta_finger_count)
species_df["bf_bool"] = species_df["bf_count"].apply(count2bool)
species_df["tree_node"] = species_df.apply(lambda x: f"{x.species} - {x.taxid}", axis=1)

tree2beta_finger = dict(zip(species_df.tree_node, species_df.bf_bool))
tree2beta_finger_count = dict(zip(species_df.tree_node, species_df.bf_count))
tree2beta_finger_count = {k: v for k,v in tree2beta_finger_count.items() if v > 0}

with open(output_beta_finger_presence, "w+") as output:
    output.write("DATASET_BINARY\n")
    output.write("SEPARATOR COMMA\n")
    output.write("DATASET_LABEL,Beta-finger Presence\n")
    output.write("FIELD_SHAPES,1\n")
    output.write("FIELD_LABELS,f1\n")

    output.write("DATA\n")
    
    for k,v in tree2beta_finger.items():
    	output.write(f"{k},{v}\n")

with open(output_poaceae_label, "w+") as output:
    output.write("DATASET_RANGE\n")
    output.write("SEPARATOR COMMA\n")
    output.write("DATASET_LABEL,Poaceae\n")
    output.write("FIELD_SHAPES,1\n")
    output.write("FIELD_LABELS,f1\n")
    output.write("LABEL_OUTLINE_WIDTH,1\n")
    output.write("LABEL_OUTLINE_COLOR,#000000\n")
    output.write("RANGE_COVER,clade\n")
    output.write("COVER_DATASETS,1\n")
    output.write("COVER_LABELS,1\n")

    output.write("DATA\n")
    output.write("Aegilops bicornis - 4483,Streptochaeta angustifolia - 38733,#ffe4c4")


text_info = "-1,#000000,italic,2,0"

with open(output_species_labels, "w+") as output:
    output.write("DATASET_TEXT\n")
    output.write("SEPARATOR COMMA\n")
    output.write("DATASET_LABEL,Monocots\n")

    output.write("DATA\n")
    
    for k,v in tree2beta_finger.items():
    	output.write(f"{k},{k.split('-')[0][:-1]},{text_info}\n")








