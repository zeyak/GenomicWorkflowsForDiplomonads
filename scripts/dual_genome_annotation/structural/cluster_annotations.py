from snakemake.shell import shell

prodigal_aa = snakemake.input.prodigal_aa
seq_identity = snakemake.params.seq_identity
length_cutoff = snakemake.params.length_cutoff

prodigal_clustered = snakemake.output.prodigal_clustered

# -c	sequence identity threshold, default 0.9
# -s	length difference cutoff, default 0.0
shell(f"cd-hit -i {prodigal_aa} -o {prodigal_clustered} -c {seq_identity} -s {length_cutoff} -T 30")
