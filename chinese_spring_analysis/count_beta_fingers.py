import pandas as pd

beta_fingers = "/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/T2T_chinese_spring_analysis/analysis/KFP_annotation/T2T_beta_finger_domains.txt"
output_annotation = "/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/T2T_chinese_spring_analysis/analysis/KFP_phylogenetics/T2T_KFP_beta_finger_counts.itol"
output_table = "/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/T2T_chinese_spring_analysis/analysis/KFP_phylogenetics/T2T_KFP_beta_finger_counts.tsv"

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

beta_finger_df = load_annotations(beta_fingers)

beta_finger_df['length'] = beta_finger_df.env_to - beta_finger_df.env_from

beta_finger_df = beta_finger_df[beta_finger_df.length >= 26]

counts = beta_finger_df.groupby(by='target').size().to_dict()

colour_per_count = {1: "#FF0000", 2: "#FFFF00", 3: "#0033CC"}

# Itol annotation of KFP beta-fingers

with open(output_annotation, "w+") as out:

    out.write("DATASET_COLORSTRIP\n")
    out.write("SEPARATOR COMMA\n")
    out.write("DATASET_LABEL,KFP_counts\n")
    out.write("HEIGHT_FACTOR,2.2\n")
    out.write("FIELD_COLORS,#FF0000,#FFFF00,#0033CC\n")
    out.write("FIELD_LABELS,1,2,3\n")
    out.write("BORDER_WIDTH,0.5\n")
    out.write("BORDER_COLOR,#000000\n")
    out.write("COMPLETE_BORDER,1\n")
    out.write("DATA\n")

    for k, v in counts.items():
        out.write(f"{k},{colour_per_count[v]},{v}\n")

# KFP counts supplementary

with open(output_table, "w+") as out:

    out.write("KFP\tbeta_finger_count\n")

    for k, v in counts.items():
        out.write(f"{k}\t{v}\n")
