#!/usr/bin/python3

KFP_positions = "/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/T2T_chinese_spring_analysis/analysis/KFP_physical_positions/KFP_loci_positions.tsv"

cloned_kfps = "/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/T2T_chinese_spring_analysis/data/KFP_genbank_ids.csv"

output_annotation = "/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/T2T_chinese_spring_analysis/analysis/KFP_phylogenetics/T2T_KFP_subgenome.itol"

# Identify subgenome for each T2T KFPs
with open(KFP_positions) as file:

	subgenome_dict = dict()

	for line in file:
		data = line.split("\t")
		subgenome = data[2][1]
		
		if subgenome not in ["A", "B", "D"]:
			print(subgenome)
		
		subgenome_dict[data[0]] = subgenome

# Generate data for cloned kfps

with open(cloned_kfps) as file:
	cloned_colour = {line.split(",")[1].strip() : "#000000" for line in file}
		

print(cloned_colour)

subgenome_colour = {'A': "#7BB662", 'B': "#FBB149", 'D': "#994FB2"}

with open(output_annotation, "w+") as out:
	
	out.write("DATASET_COLORSTRIP\n")
	out.write("SEPARATOR COMMA\n")
	out.write("DATASET_LABEL,KFP_subgenome\n")
	out.write("HEIGHT_FACTOR,2.2\n")
	out.write("FIELD_COLORS,#FF0000,#FFFF00,#0033CC\n")
	out.write("FIELD_LABELS,1,2,3\n")
	out.write("BORDER_WIDTH,0.5\n")
	out.write("BORDER_COLOR,#000000\n")
	out.write("COMPLETE_BORDER,1\n")
	out.write("DATA\n")
	
	for k, v in subgenome_dict.items():
		out.write(f"{k},{subgenome_colour[v]},{v}\n")
	
	for k, v in cloned_colour.items():
		out.write(f"{k},{v},cloned kfp\n")
	