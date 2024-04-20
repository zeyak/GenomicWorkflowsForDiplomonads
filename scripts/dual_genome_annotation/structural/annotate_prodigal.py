from snakemake.shell import shell

assembly = snakemake.input.assembly

prodigal_gff = snakemake.output.prodigal_gff
prodigal_aa = snakemake.output.prodigal_aa

shell(f"prodigal -i {assembly} -o {prodigal_gff} -g 6 -f gff -c -a {prodigal_aa}""")
