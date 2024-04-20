from snakemake.shell import shell

shell(f"diamond makedb --in diplomonad_genomes.faa --db diplomonad_genomes.db")
