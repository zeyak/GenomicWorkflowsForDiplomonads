from snakemake.shell import shell

# input
genome = snakemake.input.genome
read = snakemake.input.read
merylDB = snakemake.input.merylDB
repetitive_k15 = snakemake.input.repetitive_k15
# output
bam = snakemake.output.bam
bai = snakemake.output.bai
# params
num_threads = snakemake.params.num_threads

shell(f"""winnowmap -W {repetitive_k15} -ax map-ont {genome} {read} | samtools sort -o {bam}""")
shell(f"""samtools index {bam} {bai} """)
