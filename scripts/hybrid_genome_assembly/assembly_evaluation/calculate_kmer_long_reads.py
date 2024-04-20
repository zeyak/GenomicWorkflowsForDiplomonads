from snakemake.shell import shell

#input
genome = snakemake.input.genome

#output
merylDB= snakemake.output.merylDB
repetitive_k15 = snakemake.output.repetitive_k15

shell(f"""meryl count k=15 output {merylDB} {genome}""")
shell(f"""meryl print greater-than distinct=0.9998 {merylDB} > {repetitive_k15}""")
