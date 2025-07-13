#!/usr/bin/python3

from Bio import SeqIO
from collections import namedtuple

T2T_CS_proteome = "/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/beta_finger_HMM/data/T2T_chinese_spring/ncbi_dataset/data/GCA_040256815.2/protein.faa"

proteome = [record for record in SeqIO.parse(T2T_CS_proteome, "fasta")]

Locus = namedtuple("Locus", ["name", "length"])

loci = {}

for protein in proteome:
	
	protein_length = len(protein.seq)
	protein_name = str(protein.id)
	
	protein_locus = protein.description.split()[3]
	
	if protein_locus in loci:
		if loci[protein_locus].length < protein_length:
			loci[protein_locus] = Locus(protein_name, protein_length)
	
	else:
		loci[protein_locus] = Locus(protein_name, protein_length)


longest_protein_per_locus = {locus.name for locus in loci.values()}

longest_locus_proteome = [record for record in proteome if record.id in longest_protein_per_locus]

print(len(proteome))
print(len(longest_locus_proteome))

SeqIO.write(longest_locus_proteome, 
	"/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/T2T_chinese_spring_analysis/analysis/locus_representatives/T2T_CS_proteome.longest_gene_models.faa",
 	"fasta")

	