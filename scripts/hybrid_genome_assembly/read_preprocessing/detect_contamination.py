from snakemake.shell import shell

assembly = snakemake.input.assembly

db_prefix_human = snakemake.params.db_prefix_human
db_prefix_prok = snakemake.params.db_prefix_prok
e_value = snakemake.params.e_value
outfmt= snakemake.params.outfmt
threads = snakemake.params.threads

putative_contamination = snakemake.output.putative_contamination

#Search contamination in human database
if human == True:
    shell(f"blastn \
     -query {assembly} \
     -db {db_prefix_human} \
     -out {contamination} \
     -evalue {e_value} \
     -outfmt {outfmt} \
     -num_threads {threads}")

#Search contamination in refrence prokaryotes
else:
    shell(f"blastn \
     -query {assembly} \
     -db {db_prefix_prok} \
     -out {contamination} \
     -evalue {e_value} \
     -outfmt {outfmt} \
     -num_threads {threads}")
