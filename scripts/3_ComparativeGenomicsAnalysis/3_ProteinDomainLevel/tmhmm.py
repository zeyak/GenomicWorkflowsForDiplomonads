from snakemake.shell import shell

protein = snakemake.input.protein
out = snakemake.output[0]

shell(f"""tmhmm -f {protein} -p """)
