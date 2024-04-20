import pandas as pd

try:
    og_file = snakemake.input.og
    out = snakemake.output[0]
except NameError:
    # testing
    og_count_file = "output/1_orthofinder/Results_Oct17_2/Orthogroups/Orthogroups.GeneCount.tsv"

    og_count = "data/orthogroups/og_count.csv"

df = pd.read_csv(og_count_file, sep="\t", header='infer').rename(columns={"Orthogroup": "OG"})
df = df[["OG", "carpe", "kbiala", "HIN", "trepo", "spiro", "wb", "muris", "Total"]]

df.to_csv(og_count, sep="\t", index=False)
