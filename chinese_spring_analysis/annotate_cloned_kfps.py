#!/usr/bin/python3

cloned_kfps = "/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/T2T_chinese_spring_analysis/data/KFP_genbank_ids.csv"

output_annotation = "/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/T2T_chinese_spring_analysis/analysis/KFP_phylogenetics/cloned_kfps.itol"


# genbank2name

with open(cloned_kfps) as file:
	genbank2name = {line.split(",")[1].strip() : line.split(",")[0] for line in file}
		

with open(output_annotation, "w+") as out:
	
    out.write('DATASET_TEXT\n')
    out.write('SEPARATOR COMMA\n')
    out.write('DATASET_LABEL,Cloned KFPs\n')
    out.write("COLOR,#ff0000\n")
    out.write("MARGIN,0\n\n\n")
    out.write("DATA\n")
    
    for k, v in genbank2name.items():
    
    	if k == "Genebank_ID":
    		continue
    		
    	out.write(f"{k},{v},-1,#000000,italic,10,0\n")
	