#!/usr/bin/env python3

from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

args = sys.argv[1:]

bed = args[0] # coloured bed of kmerGWAS output
lengths_file = args[1] # can be derived from IGV index of reference genome
output = args[2] # OUTDIR 

figure = plt.figure(figsize=(12, 7),frameon=True)
ax = figure.add_subplot(111)
ax.xaxis.set_ticks_position("none")
ax.yaxis.set_ticks_position("left")
ax.set_ylabel(r'$-{log}_{10}(p)$',weight="bold")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["bottom"].set_visible(False)

lengths_dict = {}
ticks = []
annots = []

with open(lengths_file,'r') as f:
	count=0
	xmin = 0
	for l in f:
		count += 1
		lvals = l.strip().split()
		chrom = lvals[0]
		length = int(lvals[1])
		lengths_dict[chrom] = xmin
		xmax=xmin+length
		if (count+1) % 2 == 0:
			ax.axvspan(xmin=xmin, xmax=xmax, color="#E5E5E5")
		else:
			ax.axvspan(xmin=xmin, xmax=xmax, color="#c5c5c5")
		ticks.append((xmin + xmax) / 2)
		annots.append(chrom)
		xmin=xmax
	tot_len = xmax


df = pd.read_csv(bed, sep='\t', header=None)
df.columns = ['chrom', 'chromStart', 'chromEnd', 'kmer', 'p_value', 'strand', 'chromStartChunk', 'chromEndChunk','colours']

df = df[df.chrom.isin(lengths_dict.keys())]
df['corr_p_value'] = df.p_value.apply(lambda x: -np.log10(x))


def correct_positions(row) -> int:
	kmer_midpoint = (row.chromStart +row.chromEnd)/2
	return lengths_dict[row.chrom] + kmer_midpoint

df['corr_position'] = df.apply(correct_positions, axis=1)

labels = [ann[-2:] for ann in annots]

ax.set_xticks(ticks)
ax.set_xticklabels(labels)
ax.set_xlim(0,tot_len)
ax.set_ylim(df.corr_p_value.min(), df.corr_p_value.max()+1)
ax.scatter(df.corr_position, df.corr_p_value, s=8, color='#000000')
plt.savefig(output)
