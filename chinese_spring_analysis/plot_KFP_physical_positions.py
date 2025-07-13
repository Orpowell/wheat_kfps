#!/usr/bin/python3

import json

import pandas as pd 
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.transforms as mtransforms

from Bio import SeqIO
from matplotlib.patches import FancyBboxPatch
from matplotlib.patches import Patch

plt.rcParams.update({'font.size': 14})

outdir = "/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/T2T_chinese_spring_analysis/analysis/KFP_physical_positions"
sequence_report = "/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/beta_finger_HMM/data/T2T_chinese_spring/ncbi_dataset/data/GCA_040256815.2/sequence_report.jsonl"
T2T_gff = "/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/beta_finger_HMM/data/T2T_chinese_spring/ncbi_dataset/data/GCA_040256815.2/genomic.gff"
T2T_KFPs = "/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/T2T_chinese_spring_analysis/analysis/KFP_annotation/T2T_KFP_full_length.faa"
T2T_proteome = "/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/beta_finger_HMM/data/T2T_chinese_spring/ncbi_dataset/data/GCA_040256815.2/protein.faa"

# Load assembly metadata
assembly_info = pd.read_json(sequence_report, lines=True)
genbank2chr = dict(zip(assembly_info.genbankAccession, assembly_info.chrName))
chr_lengths = dict(zip(assembly_info.chrName, assembly_info.length))

# Load KFP Ids
T2T_kfps = {record.id for record in SeqIO.parse(T2T_KFPs, "fasta") if str(record.id).endswith(".1")}

locus2protein = {record.description.split()[3] : record.id for record in SeqIO.parse(T2T_proteome, 'fasta') if record.id in T2T_kfps}

# Load and preprocess GFF

gff = pd.read_csv(T2T_gff, sep="\t", skiprows=7, header=None) 
gff[0] = gff[0].replace(genbank2chr) # Convert genbank to chromosome
gff = gff[gff[2] == "gene"] # filter for genes
gff['locus'] = gff[8].apply(lambda x: x.split(';')[1][5:])
gff['protein'] = gff.locus.replace(locus2protein) # Add protein at each locus
gff = gff[gff['protein'].isin(T2T_kfps)] # Filter for KFP loci
kfp_loci = gff[['protein','locus', 0, 3, 4]]
kfp_loci[3] = kfp_loci[3].astype(int)
kfp_loci[4] = kfp_loci[4].astype(int)
kfp_loci.sort_values(by=0, inplace=True, ascending=False)
kfp_loci['length'] = kfp_loci[4] - kfp_loci[3]
kfp_loci['plot_data'] = list(zip(kfp_loci[3], kfp_loci.length))
kfp_loci.to_csv(f"{outdir}/KFP_loci_positions.tsv", sep='\t', header=None, index=False, columns=['protein','locus', 0, 3, 4])

# Count per chromosome
chr_counts = kfp_loci.groupby(by=[0]).protein.count()
chr_counts = chr_counts.to_frame().reset_index()
chr_counts.sort_values(by='protein', inplace=True, ascending=False)
chr_counts.to_csv(f"{outdir}/KFP_chromosome_counts.tsv", sep='\t', header=None, index=False)


# Count per sub-genome
chr_counts['subgenome'] = chr_counts[0].str[1]
subgenome_counts = chr_counts.groupby(by="subgenome").protein.sum()
subgenome_counts = subgenome_counts.to_frame().reset_index()
subgenome_counts.to_csv(f"{outdir}/KFP_subgenome_counts.tsv", sep='\t', header=None, index=False)


# Plot Physical Positons

ylabels = [f"Chr{chrom}"  for chrom in list(chr_lengths.keys())]
y = [n for n, chrom in enumerate(list(chr_lengths.keys()))]
x = list(chr_lengths.values())

xlabels = [i for i in range(0, 1000, 100)]

fig, ax = plt.subplots(figsize=(10, 5))

# Plot Chromosomes
for i in y:
	ax.broken_barh(xranges=[(0, x[i])], yrange=(i-0.375, 0.75), color="#bbbbbb")

# Plot Centromere for each chromosome

centromere_data = pd.read_csv('/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/T2T_chinese_spring_analysis/data/T2T_chinese_spring_centromere_coordinates.csv')

for n, data in centromere_data.iterrows():
	ax.broken_barh(xranges=[(data[1], data[2])], yrange=(n-0.375, 0.75), color="purple")

# Plot KFP loci for each chromosome

for n, data in enumerate(kfp_loci.groupby(by=[0])):
	ax.broken_barh(xranges=data[1].plot_data, yrange=(n-0.375, 0.75), color="red")


"""
# Figure Legend
legend_elements = [Patch(facecolor='red', edgecolor='black', label='KFP loci'),
					Patch(facecolor='purple', edgecolor='black', label='Centromere')]

ax.legend(handles=legend_elements, loc=0)
"""

# Axis settings

ax.set_xlim(0, 9e8)
ax.set_ylim(-1, 21)
ax.tick_params(left=False, labeltop=True, top=True, labelbottom=False, bottom=False, direction="in")
ax.set_yticks(range(21),labels=ylabels)
ax.set_xticklabels(xlabels)
ax.spines[['right', 'bottom', 'left']].set_visible(False)
ax.invert_yaxis()

ax.set_title("Mbp")

# Save
plt.tight_layout()
plt.savefig(f"{outdir}/KFP_loci_physical_positions.svg", dpi=300)
plt.show()




