from snakemake.shell import shell

glimmerhmm_annot = snakemake.input.glimmerhmm_annot
prodigal_annot = snakemake.input.prodigal_annot

merged_annot = snakemake.output.merged_annot

shell(f"MergeAnnotations.py -i {glimmerhmm_annot} {prodigal_annot} -o {merged_annot}")