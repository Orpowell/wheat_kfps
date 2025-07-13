#! /usr/bin/python3

from collections import namedtuple

from Bio import SeqIO



KFP_beta_finger_data = "/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/beta_finger_HMM/analysis/KFP_structure_analysis/KFP_beta_finger_coordinates.csv"
KFP_sequence_dir = "/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/beta_finger_HMM/data/KFP_genbank_sequences"
KFP = namedtuple("KFP", "KFP, genbank, bf_start, bf_end")



cache = {}
sequences = []

with open(KFP_beta_finger_data, "r") as file:
	
	next(file) # skip header
	
	for line in file:
		
		data = KFP(*line.strip().split(","))
		
		for record in SeqIO.parse(f"{KFP_sequence_dir}/{data.genbank}.faa", "fasta"):
			
			if record.id not in cache.keys():
				cache[record.id] = 1
			
			fasta_format = f">{record.id}_bf_{cache[record.id]}\n{str(record.seq[int(data.bf_start)-1 : int(data.bf_end)])}\n"
			
			cache[record.id] += 1
			sequences.append(fasta_format)

with open("/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/beta_finger_HMM/analysis/beta_finger_seed_HMM/KFP_beta_finger_sequences.faa", "w+") as file:
	for sequence in sequences:
		file.write(sequence)
			
			
		
		

