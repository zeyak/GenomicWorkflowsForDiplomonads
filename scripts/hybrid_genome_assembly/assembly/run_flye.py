from snakemake.shell import shell

input = snakemake.input
output = snakemake.output

# only nanopore reads
shell(f"""flye --nano-raw {input} --genome-size 114m --threads 25 -o {output}""")
