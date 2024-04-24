from snakemake.shell import shell

protein = snakemake.input.protein
out = snakemake.output

shell(f"""signalp6 --fastafile {protein} --organism eukarya --output_dir {out} --format txt,png --mode fast""")
