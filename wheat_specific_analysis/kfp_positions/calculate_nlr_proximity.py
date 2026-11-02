import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
from scipy import stats
from Bio import SeqIO

class KFP:
    def __init__(self, id, chrom, start, stop):
        self.id = id
        self.chrom = chrom
        self.start = int(start)
        self.stop = int(stop)
        self.midpoint = (self.start+self.stop) // 2

    def find_nearest_nlr(self, nlr_data):

        possible = nlr_data[(nlr_data[0] == self.chrom)]
        possible['proximity'] = abs(self.midpoint - possible.midpoint)
        possible.sort_values(by='proximity', ascending=True, inplace=True)
        self.nearest_nlr = possible.iloc[0]['proximity']


def load_kfp_data(assembly_data, gff, kfpome):

    jsonObj = pd.read_json(path_or_buf=assembly_data, lines=True)
    chr2name = dict(zip(jsonObj.genbankAccession, jsonObj.chrName))

    records = SeqIO.parse(kfpome, "fasta")

    kfp_ids = {record.id for record in records}

    KFP_data = []

    with open(gff) as file:

        for line in file:

            if line.startswith('#'):
                continue

            data = line.strip().split('\t')

            chrom = data[0]
            annotation_type = data[2]
            start = data[3]
            stop = data[4]

            if annotation_type == 'gene':

                gene_id = data[8].split('ID=')[1]

                if f"{gene_id}.1" in kfp_ids:

                    new_kfp = KFP(gene_id, chr2name[chrom], start, stop)

                    KFP_data.append(new_kfp)

    return KFP_data

def load_nlr_data(nlrs, assembly_data):

    jsonObj = pd.read_json(path_or_buf=assembly_data, lines=True)
    chr2name = dict(zip(jsonObj.genbankAccession, jsonObj.chrName))

    NLRs = pd.read_csv(nlrs, sep='\t', header=None)

    NLRs = NLRs[NLRs[2].isin(["CC-NBARC-LRR"])]

    NLRs[0] = NLRs[0].replace(chr2name)

    NLRs['midpoint'] = (NLRs[3] + NLRs[4]) // 2

    return NLRs


if __name__ == '__main__':
    
    assembly_data='GCA_040256815.2/sequence_report.jsonl' # from NCBI
    gff="wheat_chinese_spring.gff" # from Helixer
    kfpome="wheat_chinese_spring_KFP_full_length.faa" # From KFP-annotator
    nlrome="T2T_CS_NLRs.txt" # From NLR-annotator
    OUTDIR='outdir'
    
    KFP_data = load_kfp_data(assembly_data=assembly_data, gff=gff, kfpome=kfpome)

    NLR_data = load_nlr_data(nlrs=nlrome, assembly_data=assembly_data)
    
    [kfp.find_nearest_nlr(NLR_data) for kfp in KFP_data]

    nlr_proximity = [kfp.nearest_nlr for kfp in KFP_data]
    
    fig, axs = plt.subplots()
    axs.spines[["right", "top"]].set_visible(False)
    axs.set_xlabel("Distance to Nearest CC-NLR (bp)")
    axs.set_xlim(3,9)
    axs.set_ylim(0,1)
    axs.set_ylabel("Proportion of KFPs")
    sns.ecdfplot(np.log10(nlr_proximity), color='black')

    xlabels = [3,4,5,6,7,8,9]
    nu_xlabels = [f"10$^{int(x)}$" for x in xlabels]
    axs.set_xticks(xlabels, nu_xlabels)
    plt.savefig(f"{OUTDIR}/wheat_chinese_spring_nlr_proximity.svg")

    print(f"{stats.percentileofscore(a=nlr_proximity, score=150000)=}")
    print(f"{stats.percentileofscore(a=nlr_proximity, score=200000)=}")
    print(f"{stats.percentileofscore(a=nlr_proximity, score=500000)=}")

    with open(f"{OUTDIR}/wheat_chinese_spring_nlr_proximity.tsv", "w+") as file:
        file.write("KFP\tnearest CC-NLR (bp)\n")
        for kfp in KFP_data:
            file.write(f"{kfp.id}\t{kfp.nearest_nlr}\n")

