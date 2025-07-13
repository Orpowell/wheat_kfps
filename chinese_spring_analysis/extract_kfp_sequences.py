import pandas as pd
from Bio import SeqIO

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

beta_fingers = "/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/T2T_chinese_spring_analysis/analysis/KFP_annotation/T2T_beta_finger_domains.txt"
kinases = "/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/T2T_chinese_spring_analysis/analysis/KFP_annotation/T2T_pkinase_domains.txt"
sequences = "/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/T2T_chinese_spring_analysis/analysis/locus_representatives/T2T_gene_models_and_cloned_kfps.faa"
outdir = "/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/T2T_chinese_spring_analysis/analysis/KFP_annotation"

bf_df = load_annotations(beta_fingers)
bf_df['length'] = bf_df.env_to - bf_df.env_from
bf_df = bf_df[bf_df.length >= 26]

k_df = load_annotations(kinases)

merged_df = k_df.merge(bf_df, on='target')

merged_df['kinase_coords'] = pd.IntervalIndex.from_arrays(merged_df.env_from_x, merged_df.env_to_x)
merged_df['bf_coords'] = pd.IntervalIndex.from_arrays(merged_df.env_from_y, merged_df.env_to_y)
merged_df['overlap'] = [data.kinase_coords.overlaps(data.bf_coords) for index, data in merged_df.iterrows()]
merged_df = merged_df[merged_df.overlap == True]

print(merged_df[merged_df.target == 'WEM02089'])

primary_from = merged_df.groupby(by='target').env_from_x.min().to_frame()
primary_to = merged_df.groupby(by='target').env_to_x.max().to_frame()

primary_modules = primary_from.merge(primary_to, on='target')
primary_modules['coords'] = list(zip(primary_modules.env_from_x, primary_modules.env_to_x))

primary_modules_dict = primary_modules.coords.to_dict()


# write coords to tsv
with open(f'{outdir}/T2T_KFP_primary_modules.tsv', 'w+') as out:
    for k, v in primary_modules_dict.items():
        line = f"{k}\t{v[0]}\t{v[1]}\n"
        out.write(line)

# extract primary module sequences 

records = SeqIO.parse(sequences, "fasta")

seq_dict = {record.id : record.seq for record in records}

with open(f'{outdir}/T2T_KFP_primary_modules.faa', 'w+') as out:

    for k, v in primary_modules_dict.items():
        out.write(f">{k}\n")
        out.write(f"{seq_dict[k][v[0]:v[1]]}\n")

with open(f'{outdir}/T2T_KFP_full_length.faa', 'w+') as out:
    for k, v in primary_modules_dict.items():
        out.write(f">{k}\n")
        out.write(f"{seq_dict[k]}\n")