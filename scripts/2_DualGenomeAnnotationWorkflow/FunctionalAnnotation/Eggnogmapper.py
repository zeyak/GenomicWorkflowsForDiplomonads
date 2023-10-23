from snakemake.shell import shell

proteome = snakemake.input.proteome
threads = snakemake.params.threads
outdir = snakemake.params.outdir
datadir = snakemake.params.datadir

shell(f"python emapper.py -i {proteome} --output {outdir} -m diamond --data_dir {datadir}")
