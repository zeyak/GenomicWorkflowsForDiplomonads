from snakemake.shell import shell

muris = snakemake.input.muris
wb = snakemake.input.wb
spiro = snakemake.input.spiro

diplomonad_proteomes = snakemake.output.diplomonad_proteomes

shell(f"cat {muris} {wb} {spiro} > {diplomonad_proteomes.faa}")
