#! /usr/bin/python3

import json   
import glob
from Bio import SeqIO


submission_defaults = []

for file in glob.glob("/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/beta_finger_HMM/data/KFP_genbank_sequences/*.faa"):
	
	for record in SeqIO.parse(file, 'fasta'):
		name = str(record.id)
		seq = str(record.seq)
		submission = {"name": name, "modelSeeds": [], "sequences": [{"proteinChain": {"sequence": seq, "count": 1}}], "dialect": "alphafoldserver", "useStructureTemplate": True, "version": 1}
		submission_defaults.append(submission)
		
with open(f'/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/beta_finger_HMM/data/AlphaFold3_batch.json', "w+", encoding='utf-8') as file:
    json.dump(submission_defaults, file, ensure_ascii=False, indent=4)