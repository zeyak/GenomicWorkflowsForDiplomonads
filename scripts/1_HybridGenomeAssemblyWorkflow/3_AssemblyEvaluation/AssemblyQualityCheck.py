from snakemake.shell import shell

assembly = snakemake.input.assembly
outdir = snakemake.output.outdir

shell(f"quast.py {assembly} -o {outdir} --threads 2 --eukaryote")
