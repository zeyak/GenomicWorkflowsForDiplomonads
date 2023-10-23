from snakemake.shell import shell

assembly = snakemake.input.assembly
illumina_run1_R1 = snakemake.input.illumina_run1_R1
illumina_run1_R2 = snakemake.illumina_run1_R2
illumina_run2_R1 = snakemake.input.illumina_run2_R1
illumina_run2_R2 = snakemake.illumina_run2_R2
illumina_run3_R1 = snakemake.input.illumina_run3_R1
illumina_run3_R2 = snakemake.illumina_run3_R2

polished_assembly = snakemake.ouptut.polished_assembly

shell(f"""polca.sh -a {assembly} -r \
'{illumina_run1_R1} {illumina_run1_R2}
{illumina_run2_R1} {illumina_run2_R2}
{illumina_run3_R1} {illumina_run3_R2}' \
-t 30""")
