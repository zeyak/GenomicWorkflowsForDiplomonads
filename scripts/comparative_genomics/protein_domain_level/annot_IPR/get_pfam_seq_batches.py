import os
import pandas as pd
from Bio import SeqIO

# Initialize variables
pfam = "PF12799"
batch_number = 1

try:
    # Attempt to use Snakemake variables if available
    pfam = snakemake.wildcards.pfam
    gene_id_list = snakemake.input[0]
    proteome_fasta = snakemake.input[1]
    protein_fasta = snakemake.output[0]  # Assuming this isn't used since we're creating multiple output files.
except NameError:
    # Fallback for testing outside Snakemake
    gene_id_list = f"output/interproscan/processed_data/{pfam}_geneid.csv"
    proteome_fasta = "resources/Genome/HIN.faa"
    output_dir = f"output/interproscan/processed_data/{pfam}"

# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# Load the gene IDs
ids = pd.read_csv(gene_id_list, header=None, sep='\t')[0].drop_duplicates().tolist()

# Initialize a list to hold the current batch of records
current_batch_records = []

# Process the FASTA file
for record in SeqIO.parse(proteome_fasta, "fasta"):
    if record.id in ids:
        current_batch_records.append(record)

        # Check if the current batch has reached 500 records
        if len(current_batch_records) == 500:
            # Define the output filename for the current batch
            output_filename = f"{output_dir}/batch_{batch_number}.fasta"
            # Write the current batch to a new FASTA file
            with open(output_filename, "w") as output_fasta:
                SeqIO.write(current_batch_records, output_fasta, "fasta")
            # Reset the current batch records list for the next batch
            current_batch_records = []
            # Increment the batch number for the next output file
            batch_number += 1

# After the loop, check if there are any remaining records to write
if current_batch_records:
    # Define the output filename for the remaining records
    output_filename = f"{output_dir}/batch_{batch_number}.fasta"
    with open(output_filename, "w") as output_fasta:
        SeqIO.write(current_batch_records, output_fasta, "fasta")

