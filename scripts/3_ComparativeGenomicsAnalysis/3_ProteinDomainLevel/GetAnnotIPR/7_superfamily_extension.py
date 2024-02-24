import pandas as pd
import glob
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

try:
    og_file = snakemake.input.og
    out = snakemake.output[0]
except NameError:
    # testing
    ann_ipr_file = "output/interproscan/processed_data/ann_ipr_cat.csv"

    out_file1 = "output/interproscan/plot/superfamily_heatmap.png"
    out_file2 = "output/interproscan/plot/superfamily_heatmap_norm.png"
    out_file3 = "output/interproscan/processed_data/ann_ipr_cat_superfamily.csv"

def plot_heatmap(df, out_file):
    df = df.rename(columns={"carpe": "C. membranifera",
                            "kbiala": "K. bialata",
                            "HIN": "H. inflata",
                           # "trepo": "Trepomonas pc1",
                            "spiro": "S. salmonicida",
                            "wb": "G. intestinalis",
                            "muris": "G. muris"})

    df = df[["C. membranifera", "K. bialata", "H. inflata", "S. salmonicida", "G. intestinalis",
             "G. muris"]]
    custom_labels = ["LRR", "HB", "CP", "PK", "CRP", "WD40", "ankyrin"]

    sns.set_theme(style="whitegrid")

    f, ax = plt.subplots(figsize=(10, 10))

    heatmap = sns.heatmap(df.T,
                     norm=LogNorm(),
                     cmap=sns.color_palette("ch:start=.2,rot=-.3", as_cmap=True),
                     square=True,
                     fmt='g',
                     linewidths=.4,
                     annot=True,
                     cbar_kws={"shrink": .5}
                     )
    #heatmap.set_title("top7_superfamily_normalized")
    heatmap.set_xlabel('')  # Removes the x-axis title

    # Positioning custom labels above each column
    for idx, label in enumerate(custom_labels):
        ax.text(x=idx + 0.5, y=-0.1, s=label, ha='center', va='center')

    plt.savefig(out_file, format="png", bbox_inches='tight', dpi=1200)
    plt.show()
    return plt

def dict_count_ipr(df, sp) -> pd.DataFrame:
    """
    Count iprs for superfamily in each sp
    @return: superfamily count table, dataframe
    """
    df_ = df.groupby(["sp"]).get_group(sp)
    df_ = df_[df_["ann_inter"].str.contains("superfamily", case=False)]
    count = df_.groupby(["ipr", "ann_inter"]).size().reset_index().sort_values(by=[0], ascending=False)
    return count.rename(columns={0: sp}), df_


def normalization(df):
    df_c = df.copy()  # Make a copy of the input DataFrame
    df_c["HIN"] = (df_c["HIN"] / 142.6).round(1)
    df_c["spiro"] = (df_c["spiro"] / 14.7).round(1)
    df_c["wb"] = (df_c["wb"] / 12.6).round(1)
    df_c["muris"] = (df_c["muris"] / 9.8).round(1)
    df_c["carpe"] = (df_c["carpe"] / 24.2).round(1)
    df_c["kbiala"] = (df_c["kbiala"] / 51).round(1)
    # trepo doesn't have a genome size
    return df_c

species = ['HIN', 'spiro', 'wb', 'muris', 'carpe', 'kbiala'] #trepo is excluded

df = pd.read_csv(ann_ipr_file, header="infer", sep="\t")
df = df.dropna(subset="ipr").drop_duplicates(subset=["id", "ipr"])
count, df_superfamily = dict_count_ipr(df, "HIN")

for sp in species:
    count = pd.merge(count, dict_count_ipr(df, sp)[0], how="outer")
    df_superfamily = pd.merge(df_superfamily, dict_count_ipr(df, sp)[1], how="outer")


counts = count.set_index("ipr")
counts_norm = normalization(counts)
print(counts)

plot_heatmap(counts[:7], out_file1) # top 7 superfamily
plot_heatmap( counts_norm[:7], out_file2)
df_superfamily.to_csv(out_file3, sep="\t", index=False)