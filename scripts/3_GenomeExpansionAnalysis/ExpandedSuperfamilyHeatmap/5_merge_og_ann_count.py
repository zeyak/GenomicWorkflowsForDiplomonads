import pandas as pd

try:
    og_file = snakemake.input.og
    out = snakemake.output[0]
except NameError:
    # testing
    og_stack = "data/orthogroups/og_stack.csv"
    og_count = "data/orthogroups/og_count.csv"
    og_count_s = "data/orthogroups/og_count_s.csv"
    ann_f_cat = "data/fasta_ann/concat/ann_f_cat.csv"
    ann_ipr_cat = "data/interpro_ann/concat/ann_ipr_cat.csv"
    ann_egg_cat = "data/eggnog_ann/ann_egg_cat.csv"

    out_file = "data/orthogroups/og_ann.csv"

og_stack_ = pd.read_csv(og_stack, sep="\t", header='infer')
og_count_s_ = pd.read_csv(og_count_s, sep="\t", header='infer')  # Could be also og_count

ann_f_cat_ = pd.read_csv(ann_f_cat, sep="\t", header='infer')
ann_ipr_cat_ = pd.read_csv(ann_ipr_cat, sep="\t", header='infer')
ann_egg_cat_ = pd.read_csv(ann_egg_cat, sep="\t", header='infer')

df = pd.merge(og_stack_, og_count_s_)

df = pd.merge(df, ann_f_cat_, how="left", on="id")
df = pd.merge(df, ann_ipr_cat_, how="left")
df = pd.merge(df, ann_egg_cat_, how="left")

df = df[["OG", "id", "ann_f", "ann_inter", "ann_egg", "db", "db_acc", "ann_db", "ipr",
         "COG", "KEGG_KOs", "carpe", "kbiala", "HIN",
         "trepo", "spiro", "wb", "muris", "Total", "Type", "sp"]]





# df = get_ann(get_og_genes(og_stack, values))
# df = mark_trepo_lgt(df)

df.to_csv(out_file, sep="\t", index=False)
# df.to_excel(outfile_xlsx.format(keys), index=False)
