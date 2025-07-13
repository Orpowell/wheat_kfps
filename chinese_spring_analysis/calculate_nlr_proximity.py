#!/usr/bin/python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
from scipy import stats

pd.options.mode.copy_on_write = True

NLR_positions = "/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/T2T_chinese_spring_analysis/analysis/NLR_proximity/T2T_CS_NLRs.txt"
KFP_positions = "/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/T2T_chinese_spring_analysis/analysis/KFP_physical_positions/KFP_loci_positions.tsv"
KFP_auxiliary_modules = "/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/T2T_chinese_spring_analysis/analysis/KFP_phylogenetics/T2T_KFP_auxiliary_modules.tsv"
sequence_report = "/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/beta_finger_HMM/data/T2T_chinese_spring/ncbi_dataset/data/GCA_040256815.2/sequence_report.jsonl"

output_graph = "/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/T2T_chinese_spring_analysis/analysis/NLR_proximity/ECDF_NLR_distance_plot.svg"
output_table = "/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/T2T_chinese_spring_analysis/analysis/NLR_proximity/T2T_KFP_NLR_distance.tsv"

class KFP:
	
	def __init__(self, name, chrom, start, end, aux):
		self.name = name
		self.chrom = chrom
		self.start = int(start)
		self.end = int(end)
		self.midpoint = (start+end) // 2
		self.aux = aux
	
	def find_NLR_candidate(self, nlr_data: pd.DataFrame):
		
		possible = nlr_data[nlr_data[0] == self.chrom]
		possible['midpoint'] = (possible[3] + possible[4]) // 2		
		possible['proximity'] = abs(possible.midpoint - self.midpoint)
		possible = possible.sort_values(by='proximity', ascending=True)
		self.nearest_NLR_distance = possible.iloc[0]['proximity']
		self.nearest_NLR =  possible.iloc[0][1]
		
# Load and remove cloned KFPs

KFPs = pd.read_csv(KFP_positions, header=None, sep="\t")
aux = pd.read_csv(KFP_auxiliary_modules, header=None, sep="\t")
aux = aux[aux[0].str.endswith(".1")]
CS_KFPs = pd.merge(KFPs[[0,2,3,4]], aux, on=[0])

# Load NLRs and correct chromosome names

assembly_info = pd.read_json(sequence_report, lines=True)
genbank2chr = dict(zip(assembly_info.genbankAccession, assembly_info.chrName))
NLRs = pd.read_csv(NLR_positions, sep='\t', header=None)

NLRs[0] = NLRs[0].replace(genbank2chr)

def check_overlap(row):
	if True in KFP_intervals.overlaps(pd.Interval(int(row[3]), int(row[4]), closed="both")):
		return True
	else:
		return False

NLR_KFPs = CS_KFPs[CS_KFPs[1].isin(["NLR", "NBARC-LRR", "LRR"])]
KFP_intervals = pd.arrays.IntervalArray.from_tuples(list(zip(NLR_KFPs[3], NLR_KFPs[4])))
NLRs['KFP'] = NLRs.apply(check_overlap, axis=1)
CS_NLRs = NLRs[(NLRs[2].isin(["CC-NBARC-LRR"])) & (~NLRs.KFP)]
CS_NLRs.reset_index(drop=True, inplace=True)



# Find nearest NLR to each KFP
kfp_nlr_distance = {}
kfp_nlr = {}

for n, data in CS_KFPs.iterrows():
	kfp = KFP(name=data[0],chrom=data[2], start=data[3], end=data[4], aux=data[1])
	kfp.find_NLR_candidate(CS_NLRs)
	kfp_nlr_distance[kfp.name] = int(kfp.nearest_NLR_distance)
	kfp_nlr[kfp.name] = kfp.nearest_NLR

CS_KFPs['nlr_distance'] = CS_KFPs[0].map(kfp_nlr_distance)
CS_KFPs['nlr'] = CS_KFPs[0].map(kfp_nlr)
CS_KFPs = CS_KFPs.sort_values(by='nlr_distance', ascending=True)
CS_KFPs['cumulative'] = np.linspace(1, len(CS_KFPs), len(CS_KFPs))
CS_KFPs.to_csv(output_table, sep="\t", columns=[0, 'nlr_distance'], header=None, index=False)

# Plot 
fig, axs = plt.subplots()
axs.spines[["right", "top"]].set_visible(False)
axs.set_xlabel("Distance to Nearest CC-NLR (bp)")
axs.set_xlim(3,9)
axs.set_ylim(0,1)
axs.set_ylabel("# of KFPs")
sns.ecdfplot(np.log10(CS_KFPs['nlr_distance']), color='black')

xlabels = [3,4,5,6,7,8,9]
nu_xlabels = [f"10$^{int(x)}$" for x in xlabels]
axs.set_xticks(xlabels, nu_xlabels)
plt.savefig(output_graph)


print(f"{stats.percentileofscore(a=CS_KFPs['nlr_distance'], score=150000)=}")
print(f"{stats.percentileofscore(a=CS_KFPs['nlr_distance'], score=200000)=}")