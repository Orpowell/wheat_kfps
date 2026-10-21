import pandas as pd

monocot_jsonl_path = "/jic/scratch/groups/Paul-Nicholson/ollie/rmg2/data/monocot_origin/genomes/monocot_4447_reference_genomes/ncbi_dataset/data/assembly_data_report.jsonl"
s_angustfolia_jsonl_path = "/jic/scratch/groups/Paul-Nicholson/ollie/rmg2/data/monocot_origin/genomes/s_angustfolia/ncbi_dataset/data/assembly_data_report.jsonl"
output_table = "/jic/scratch/groups/Paul-Nicholson/ollie/rmg2/analysis/monocot_origin/monocot_genome_species.tsv"

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
