from snakemake.shell import shell

pacbio = snakemake.input.pacbio
nanopore = snakemake.input.nanopore
run1_R2 = snakemake.input.run1_R2 #take as single reads
run2 = snakemake.input.run2
run3 = snakemake.input.run3

out = snakemake.output.out
outraw = snakemake.output.outraw
num_threads = snakemake.params.num_threads

ill1_S = snakemake.params.ill1_S
ill2_P = snakemake.params.ill2_P
ill3_P = snakemake.params.ill3_P
pac = snakemake.params.pac
nano = snakemake.params.nano

shell(f"""plotCoverage -b {run1_R2} {run2} {run3} {pacbio} {nanopore}\
 --labels {ill1_S} {ill2_P} {ill3_P} {pac} {nano}\
  -o {out} -p {num_threads} --verbose --outRawCounts {outraw} """)
