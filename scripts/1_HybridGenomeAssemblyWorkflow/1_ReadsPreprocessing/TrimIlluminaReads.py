from snakemake.shell import shell

r1 = snakemake.input[0]
r2 = snakemake.input[1]
r1_u = snakemake.output[0] #unique reads
r2_u = snakemake.output[1]
r1_d = snakemake.output[2] #duplicate reads
r2_d = snakemake.output[3]
threads = snakemake.params.threads

shell(f"""trimmomatic PE -threads {threads} {r1} {r2} {r1_u} {r1_d} {r2_u} {r2_d} CROP:100 MINLEN:50""")

