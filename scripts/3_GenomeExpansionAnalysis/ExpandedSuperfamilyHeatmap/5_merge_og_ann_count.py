import pandas as pd

try:
    og_file = snakemake.input.og
    out = snakemake.output[0]
except NameError:
    # testing
    og_stack = "data/orthogroups/og_stack.csv" # Orthogroups
    og_count = "data/orthogroups/og_count.csv" # Orthogroups with gene count per species
    og_count_s = "data/orthogroups/og_count_s.csv" # Orthogroups and singletons with gene count per species
    ann_f_cat = "data/fasta_ann/concat/ann_f_cat.csv" # Annotation by NCBI fasta
    ann_ipr_cat = "data/interpro_ann/concat/ann_ipr_cat.csv" # Annotation by InterProScan
    ann_egg_cat = "data/eggnog_ann/ann_egg_cat.csv" # Annotation by eggNOG-mapper

    out_file = "data/orthogroups/og_ann.csv" # Orthogroups and singletons with gene count per species and annotation


og_stack_ = pd.read_csv(og_stack, sep="\t", header='infer')
# Use og_count if you want to exclude singletons
og_count_s_ = pd.read_csv(og_count_s, sep="\t", header='infer')

ann_f_cat_ = pd.read_csv(ann_f_cat, sep="\t", header='infer')
ann_ipr_cat_ = pd.read_csv(ann_ipr_cat, sep="\t", header='infer')
ann_egg_cat_ = pd.read_csv(ann_egg_cat, sep="\t", header='infer')

df = pd.merge(og_stack_, og_count_s_, how="inner", on="id")
df = pd.merge(df, ann_f_cat_, how="left", on="id")
df = pd.merge(df, ann_ipr_cat_, how="left", on="id")
df = pd.merge(df, ann_egg_cat_, how="left", on="id")

df = df[["OG", "id", "ann_f", "ann_inter", "ann_egg", "db", "db_acc", "ann_db", "ipr",
         "COG", "KEGG_KOs", "carpe", "kbiala", "HIN",
         "trepo", "spiro", "wb", "muris", "Total", "Type", "sp"]]


df.to_csv(out_file, sep="\t", index=False)
# df.to_excel(outfile_xlsx.format(keys), index=False)
