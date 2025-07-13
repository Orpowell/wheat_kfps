import pandas as pd

monocot_jsonl_path = "/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/monocot_KFP_presence/data/monocot_assemblies.jsonl"
s_angustfolia_jsonl_path = "/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/monocot_KFP_presence/data/s_angustfolia.jsonl"
output_table = "/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/monocot_KFP_presence/analysis/monocot_genome_species.tsv"

def load_species_data(jsonl_path: str) -> pd.DataFrame:
	
	jsonObj: pd.DataFrame = pd.read_json(path_or_buf=jsonl_path, lines=True, orient='records')

	jsonObj['taxId'] = jsonObj.assemblyInfo.apply(lambda x: x['biosample']['description']['organism']['taxId'])

	jsonObj['species'] = jsonObj.assemblyInfo.apply(lambda x: x['biosample']['description']['organism']['organismName'])

	species_data = jsonObj[['species','accession', 'taxId']]
	
	return species_data

monocots = load_species_data(monocot_jsonl_path)
s_angustfolia = load_species_data(s_angustfolia_jsonl_path)	
species_data = pd.concat([monocots, s_angustfolia])

species_data.to_csv(output_table, sep="\t", header=None, index=False)