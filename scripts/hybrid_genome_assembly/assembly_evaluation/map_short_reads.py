from snakemake.shell import shell
import os

bam = snakemake.output.bam
bai = snakemake.output.bai

num_threads = snakemake.params.num_threads
index = os.path.commonprefix(snakemake.input.index).rstrip(".")
paired = snakemake.params.get('paired', False)

if paired:
    illumina_R1 = snakemake.input.illumina_R1
    illumina_R2 = snakemake.input.illumina_R2
    shell(f"""bowtie2 -p {num_threads} -1 {illumina_R1} -2 {illumina_R2} -x {index} | samtools sort -o {bam}""")
    shell(f"""samtools index {bam} {bai}""")
else:
    read = snakemake.input.read
    shell(f"""bowtie2 -p {num_threads} {read} -x {index} | samtools sort -o {bam}""")
    shell(f"""samtools index {bam} {bai}""")