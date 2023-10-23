from snakemake.shell import shell

trepo_proteome= snakemake.input.trepo_proteome

db_prefix = snakemake.params.db_prefix
outfmt = snakemake.params.outfmt
threads = snakemake.params.threads
max_target_seqs = snakemake.params.max_target_seqs
max_hsps = snakemake.params.max_hsps
more_sensitive = snakemake.params.more_sensitive

output = snakemake.output

shell(f"""diamond blastp --query {trepo_proteome} \
 --db {db_prefix} \
 --out {output} \
 --outfmt {outfmt} \
 -threads {threads} \
 -max_target_seqs {max_target_seqs} \
 -max_hsps {max_hsps} \
 --more-sensitive {more_sensitive}
""")