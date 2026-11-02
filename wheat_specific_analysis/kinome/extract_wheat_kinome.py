import pandas as pd
from Bio import SeqIO

beta_fingers="wheat_chinese_spring_beta_finger_domains.txt" # from KFP-Annotator
kinases = "wheat_chinese_spring_pkinase_domains.txt" # from KFP-Annotator
proteome = "wheat_chinese_spring.faa" # from Helixer

OUTDIR='outdir'
output_beta_finger_presence=f'{OUTDIR}/wheat_chinese_spring_kinome_beta_fingers.itol'
output_wheat_chinese_spring_kinome=f'{OUTDIR}/wheat_chinese_spring_kinome.faa'

def load_annotations(annotations) -> pd.DataFrame:

    with open(annotations) as file:
        annot_array = []
        for line in file:
            if line.startswith("#"):
                continue

            row = [data for data in line.split(" ") if data != ''][:22]

            annot_array.append(row)

    labels = ["target", "accession", "tlen", "query", "accession", "qlen", "E-value", "score", "bias", "n", "total", "c-Evalue", "i-Evalue", "score", "bias", "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to", "acc"]

    df = pd.DataFrame(annot_array, columns=labels)

    df.env_from = df.env_from.astype(int)

    df.env_to = df.env_to.astype(int)

    return df[['target', 'query', "env_from", "env_to"]]

def load_proteome(sequences) -> dict:
    
    records = SeqIO.parse(sequences, "fasta")

    return {record.id : record.seq for record in records}


bf_df = load_annotations(beta_fingers)
bf_df['length'] = bf_df.env_to - bf_df.env_from
bf_df = bf_df[bf_df.length >= 26]

k_df = load_annotations(kinases)

merged_df = k_df.merge(bf_df, on='target')

merged_df['kinase_coords'] = pd.IntervalIndex.from_arrays(merged_df.env_from_x, merged_df.env_to_x)
merged_df['bf_coords'] = pd.IntervalIndex.from_arrays(merged_df.env_from_y, merged_df.env_to_y)
merged_df['overlap'] = [data.kinase_coords.overlaps(data.bf_coords) for index, data in merged_df.iterrows()]

bfk_df = merged_df[merged_df.overlap == True]

k_df['header'] = k_df.apply(lambda x: f"{x.target}_{x.env_from}_{x.env_to}", axis=1)
bfk_df['header'] = bfk_df.apply(lambda x: f"{x.target}_{x.env_from_x}_{x.env_to_x}", axis=1)

bfk_dict = dict(zip(bfk_df.header, bfk_df.overlap))

k_df['beta_finger_present'] = k_df['header'].map(bfk_dict)

k_df['beta_finger_present'].fillna(False, inplace=True)

bf_presence = dict(zip(k_df.header, k_df.beta_finger_present))

# Generate iTOL annotation for beta-finger containing kinase domains

with open(output_beta_finger_presence, "w+") as output:
    output.write("DATASET_COLORSTRIP\n")
    output.write("SEPARATOR COMMA\n")
    output.write("DATASET_LABEL,Beta-finger Presence\n")
    output.write("FIELD_SHAPES,2\n")
    output.write("HEIGHT_FACTOR,2.2\n")
    output.write("FIELD_LABELS,beta_finger_presence\n")
    output.write("DATA\n")
    
    for k,v in bf_presence.items():
    	if v:
    		output.write(f"{k},#dc267f,1\n")
              
# Extract kinase domains from proteome

proteome_seqs = load_proteome(proteome)

i = 0
with open(output_wheat_chinese_spring_kinome, "w+") as output:

    for n, row in k_df.iterrows():

        try:
            seq = proteome_seqs[row.target]

        except KeyError:
            print(f'{row.target} not in CS proteome')
            continue
        
        if len(seq) < 100:
            continue

        header = row.header
        domain = seq[row.env_from-1 : row.env_to-1]

        output.write(f'>{header}\n')
        output.write(f'{domain}\n')
        i+=1

print(f'{i} kinase domains written...')
