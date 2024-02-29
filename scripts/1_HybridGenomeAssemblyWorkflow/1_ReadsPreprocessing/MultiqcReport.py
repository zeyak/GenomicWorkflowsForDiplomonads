from snakemake.shell import shell

reads = snakemake.input.reads
outdir = snakemake.output.outdir

shell(f"multiqc --outdir {outdir}")