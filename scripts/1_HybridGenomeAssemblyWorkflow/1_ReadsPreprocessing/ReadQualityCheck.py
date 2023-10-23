from snakemake.shell import shell

reads = snakemake.input.reads
outdir = snakemake.output.outdir

shell(f"fastqc --threads 8 {input} --outdir {outdir}")

