from snakemake.shell import shell

protein = snakemake.input.protein
out = snakemake.output[0]

shell(f"""tmhmm -m TMHMM2.0.model -f {protein} -p """)
