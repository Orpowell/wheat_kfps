import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.transforms as mtransforms
from matplotlib.patches import FancyBboxPatch
from matplotlib.patches import Patch
from Bio import SeqIO

class KFP:
    def __init__(self, id, chrom, start, stop):
        self.id = id
        self.chrom = chrom
        self.start = start
        self.stop = stop

    def get_subgenome(self):
        return self.chrom[1]


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


def generate_itol_subgenome_annotation(KFPs, OUTDIR):

    subgenome_colour = {'A': "#CCBB44", 'B': "#4477AA", 'D': "#66CCEE"}

    with open(f"{OUTDIR}/wheat_chinese_spring_KFP_subgenomes.itol", "w+") as out:
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

        for kfp in KFPs:

            out.write(f"{kfp.id}.1,{subgenome_colour[kfp.get_subgenome()]},{kfp.get_subgenome()}\n")


def load_centromere_data(centromere_data):

    with open(centromere_data) as file:

        next(file) # skip header

        return [line.strip().split(',') for line in file]
            

def plot_kfp_positions(KFPs, assembly_data, centromere_positions):

    jsonObj = pd.read_json(path_or_buf=assembly_data, lines=True)

    chr_lengths = dict(zip(jsonObj.chrName, jsonObj.length))

    chr_n = {chrom : n for n, chrom in enumerate(jsonObj.chrName.unique())}

    print(chr_n)

    plt.rcParams.update({'font.size': 14})
    fig, ax = plt.subplots(figsize=(10, 5))

    ax.set_xlim(0, 9e8)
    ax.set_ylim(-1, 21)
    ax.tick_params(left=False, labeltop=True, top=True, labelbottom=False, bottom=False, direction="in")
    ax.set_yticks(range(21),labels=chr_lengths.keys())

    labels = [0, 100, 200, 300, 400, 500, 600, 700, 800, "     900 Mbp"]

    ax.set_xticklabels(labels)
    ax.spines[['right', 'bottom', 'left']].set_visible(False)
    ax.invert_yaxis()

    for n, length in enumerate(chr_lengths.values()):
        ax.broken_barh(xranges=[(0, length)], yrange=(n-0.375, 0.75), color="#bbbbbb")


    for centromere in centromere_positions:

        ax.broken_barh(xranges=[(int(centromere[1]), int(centromere[2]))], yrange=(chr_n[centromere[0]]-0.375, 0.75), color="purple")


    for kfp in KFPs:

        kfp_length = int(kfp.stop) - int(kfp.start)
        ax.broken_barh(xranges=[(int(kfp.start), kfp_length)], yrange=(chr_n[kfp.chrom]-0.375, 0.75), color="red")

    plt.savefig(f"{OUTDIR}/wheat_chinese_spring_KFP_positions.svg")

if __name__ == '__main__':

    assembly_data='GCA_040256815.2/sequence_report.jsonl' # from NCBI
    gff="wheat_chinese_spring.gff" # from Helixer
    kfpome="wheat_chinese_spring_KFP_full_length.faa" # from KFP-Annotator
    centromere_data='T2T_chinese_spring_centromere_coordinates.csv' # from original T2T paper supplementary
    OUTDIR='outdir'

    KFP_data = load_kfp_data(assembly_data=assembly_data, gff=gff, kfpome=kfpome)

    generate_itol_subgenome_annotation(KFPs=KFP_data, OUTDIR=OUTDIR)

    centromere_positions = load_centromere_data(centromere_data=centromere_data)

    plot_kfp_positions(KFPs=KFP_data, assembly_data=assembly_data, centromere_positions=centromere_positions)
