from snakemake.shell import shell

config = snakemake.input.config
output = snakemake.output

# MaSuRCA (FLYE_ASSEMBLY = 1, SOAP_ASSEMBLY=0)
shell(f"masurca {config}")
shell(f"bash assemble.sh > {output}")
