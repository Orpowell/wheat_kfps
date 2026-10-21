import pandas as pd

taxa_data = "/jic/scratch/groups/Paul-Nicholson/ollie/rmg2/analysis/monocot_origin/moncot_tax_id.tsv"
monocot_species = "/jic/scratch/groups/Paul-Nicholson/ollie/rmg2/analysis/monocot_origin/monocot_genome_species.tsv"
family_annotation_output = "/jic/scratch/groups/Paul-Nicholson/ollie/rmg2/analysis/monocot_origin/monocot_families.itol"
monocot_table_output = "/jic/scratch/groups/Paul-Nicholson/ollie/rmg2/analysis/monocot_origin/monocot_presence.tsv"

taxa_df = pd.read_csv(taxa_data, sep="\t")
monocot_df = pd.read_csv(monocot_species, sep="\t", header=None)

taxa_id = monocot_df[2].to_list()

species_of_interest = taxa_df[taxa_df.Taxid.isin(taxa_id)]

family_values = species_of_interest['Family name'].value_counts()

min_count = 10 

families = {k: v for k,v in family_values.items() if v >= min_count}
print(len(families))

annotated_species = species_of_interest[species_of_interest['Family name'].isin(families.keys())]

colours = [
    "#3cb44b",
    "#ffe119",
    "#4363d8",
    "#f58231",
    "#e6194b",
    "#911eb4",
    "#46f0f0",
    "#f032e6",
    "#bcf60c",
    "#fabebe",
    "#008080",
    "#e6beff",
    "#9a6324",
    "#fffac8",
    "#800000",
    "#aaffc3",
    "#ffd8b1",
    "#000075",
]

colour_dict = dict(zip(annotated_species['Family name'].unique(), colours))

species_of_interest['node_id'] = species_of_interest.apply(lambda x: f"{x['Tax name']} - {x['Taxid']}" , axis=1)

species_of_interest['family_colour'] = species_of_interest['Family name'].map(colour_dict).fillna("#ffffff")

annotation_dict = dict(zip(species_of_interest['Tax name'], species_of_interest.family_colour))

for k,v in colour_dict.items():
	print(k, v)


with open(family_annotation_output, "w+") as out:

    out.write("DATASET_COLORSTRIP\n")
    out.write("SEPARATOR COMMA\n")
    out.write("DATASET_LABEL,KFP_Annotations\n")
    out.write("HEIGHT_FACTOR,5\n")
    out.write("BORDER_WIDTH,0.5\n")
    out.write("BORDER_COLOR,#000000\n")
    out.write("COMPLETE_BORDER,1\n")
    out.write("STRIP_WIDTH,100\n")
    out.write("STRIP_LABEL_AUTO_COLOR,1\n")
    out.write("DATA\n")
    
    for k,v in annotation_dict.items():
    	out.write(f"'{k}',{v}\n")

species_of_interest.to_csv(monocot_table_output, sep="\t", columns=['Tax name', 'Taxid', 'Family name'], index=None)
