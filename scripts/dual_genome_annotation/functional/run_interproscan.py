from snakemake.shell import shell

proteome = snakemake.input.proteome
threads = snakemake.params.threads
outdir = snakemake.params.outdir

# remove the stars representing the stop codons
# sed 's/*//' trepo.faa > trepo_aa.fasta
shell(f"""interproscan.sh -i {proteome} -d {outdir} -f gff3,tsv,json -iprlookup -goterms --pathways -cpu {threads} """)
