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
output_beta_finger_table="/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/monocot_KFP_presence/analysis/monocot_beta_finger_presence.tsv"

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
    output.write("FIELD_SHAPES,2\n")
    output.write("FIELD_LABELS,beta_finger_presence\n")
    output.write("FIELD_COLORS,#000000\n")
    output.write("DATA\n")
    
    for k,v in tree2beta_finger.items():
    	if v:
    		output.write(f"{k},{v}\n")

species_df.to_csv(output_beta_finger_table, sep="\t", columns=["species", "taxid", "genome", 'bf_bool'], index=None)





