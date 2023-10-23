from snakemake.shell import shell

outdir = snakemake.params.outdir

shell(f"multiqc {outdir}")