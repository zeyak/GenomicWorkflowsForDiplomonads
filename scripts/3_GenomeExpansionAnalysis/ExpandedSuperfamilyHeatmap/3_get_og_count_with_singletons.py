import pandas as pd

try:
    og_file = snakemake.input.og
    out = snakemake.output[0]
except NameError:
    # testing
    og_count_file = "output/1_orthofinder/Results_Oct17_2/Orthogroups/Orthogroups.GeneCount.tsv"
    singleton_file = "output/1_orthofinder/Results_Oct17_2/Orthogroups/Orthogroups_UnassignedGenes.tsv"

    og_count_s = "data/orthogroups/og_count_s.csv"

df_count = pd.read_csv(og_count_file, sep="\t", header='infer').rename(columns={"Orthogroup": "OG"})
df_count = df_count[["OG", "carpe", "kbiala", "HIN", "trepo", "spiro", "wb", "muris", "Total"]]
df_count = df_count.set_index("OG").sort_values(by="Total", ascending=False)
df_count.loc[df_count["Total"] > 1, "Type"] = "OG"

# "Single genes which are excluded from the OG"
df_sing = pd.read_csv(singleton_file, sep="\t", header='infer').rename(columns={"Orthogroup": "OG"})
df_sing = df_sing.set_index("OG").fillna(0)
df_sing["Total"] = 1
df_sing = df_sing.applymap(lambda x: 1 if isinstance(x, str) == True else x)

# "Concatanate OG and singleton dataframes"
df_count_s = pd.concat([df_count, df_sing], axis=0)
df_count_s.loc[df_count_s["Total"] > 1, "Type"] = "OG"
df_count_s.loc[df_count_s["Total"] == 1, "Type"] = "singleton"

df_count_s.reset_index().to_csv(og_count_s, sep="\t", index=False)
