#! /usr/bin/python3

from Bio import SeqIO
from Bio import Entrez
from Bio import GenBank


class cloned_gene:
    def __init__(self, genbank, sequence) -> None:
        self.genbank = genbank
        self.sequence = sequence

    def __repr__(self) -> str:
        return f">{self.genbank}\n{self.sequence}\n"
    
    def save_gene_as_fasta(self, directory: str):
    	with open(f"{directory}/{self.genbank}.faa", "w+") as file:
    		file.write(self.__repr__())


if __name__ == "__main__":

	Entrez.email = "oliver.powell@kaust.edu.sa"

	input_file = "/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/beta_finger_HMM/data/KFP_genbank_ids.csv"

	with open(input_file) as file:
		genes = [line.split(",")[1].strip() for line in file]

	genbank_ids = ",".join(genes)
	
	handle = Entrez.efetch(db="protein", id=genbank_ids, rettype="gb", retmode="text")
	
	records = GenBank.parse(handle)
	
	genbank_data = [cloned_gene(genbank=record.accession[0],sequence=record.sequence) for record in records]
	
	for gene in genbank_data:
		gene.save_gene_as_fasta("/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/beta_finger_HMM/data/KFP_genbank_sequences")
		
