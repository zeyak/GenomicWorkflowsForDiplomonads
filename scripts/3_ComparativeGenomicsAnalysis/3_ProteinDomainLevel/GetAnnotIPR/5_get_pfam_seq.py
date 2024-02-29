import pandas as pd
from Bio import SeqIO

pfam= "PF12799"
try:
    pfam = snakemake.wildcards.pfam
    gene_id_list = snakemake.input[0]
    proteome_fasta = snakemake.input[1]
    protein_fasta = snakemake.output[0]
except NameError:
    # testing
    gene_id_list = f"""output/interproscan/processed_data/{pfam}_geneid.csv"""
    proteome_fasta = "resources/Genome/HIN.faa"
    protein_fasta= f"""output/interproscan/processed_data/{pfam}.faa"""

#get hit ids from blastn
id = pd.read_csv(gene_id_list, header=None, sep='\t')[0]#.tolist()
print("id", len(id))
id_list= id.drop_duplicates().tolist()

#extract fasta file for missing transcripts
with open(proteome_fasta, "r") as fasta_in:
    with open(protein_fasta, "w") as fasta_out:
        for record in SeqIO.parse(fasta_in, "fasta"):
            if record.id in id_list:
                # record.description = record.id
                SeqIO.write(record, fasta_out, "fasta")