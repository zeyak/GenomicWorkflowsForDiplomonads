import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.colors import LogNorm

"""
This script;
Gives the number of GENES with PFAM domains for the Top 7 superfamily
Makes a heatmap with
Top 10 and Top 5 PFAM domains
"""

ann_ipr_file = "output/interproscan/processed_data/ann_ipr_cat.csv"  # get pfam domains
top7sf_file = "output/interproscan/processed_data/ann_ipr_cat_top7sf.csv"  # get genes from top 7 superfamily
out_file2 = "output/interproscan/plot/superfamily_pfam_heatmap_top10.png"
out_file3 = "output/interproscan/plot/superfamily_pfam_heatmap_top5.png"
out_file4 = "output/interproscan/plot/superfamily_pfam_heatmap_top20.png"


def plot_heatmap(df, out_file):
    df = df.set_index(["db_acc", "family"])

    df = df.rename(columns={"carpe": "C. membranifera",
                            "kbiala": "K. bialata",
                            "HIN": "H. inflata",
                            "trepo": "Trepomonas pc1",
                            "spiro": "S. salmonicida",
                            "wb": "G. intestinalis",
                            "muris": "G. muris"})

    df = df[["C. membranifera", "K. bialata", "H. inflata", "Trepomonas pc1", "S. salmonicida", "G. intestinalis",
             "G. muris"]]

    sns.set_style("whitegrid", {'axes.grid': False})

    plt.figure(figsize=(5, 15))
    # Create the heatmap
    ax = sns.heatmap(df,
                     norm=LogNorm(),
                     # cmap=sns.color_palette("light:b", as_cmap=True),
                     cmap=sns.color_palette("gray_r", as_cmap=True),
                     # cmap=sns.color_palette("ch:start=.2,rot=-.3", as_cmap=True),
                     square=True,
                     fmt='g',
                     linewidths=0.3,
                     linecolor="whitesmoke",
                     annot=True,
                     cbar_kws={"shrink": .3},
                     annot_kws={"size": 6})  # Adjust font size of annotations here

    ax.tick_params(axis='both', which='major', labelsize=8)  # Adjust the size as needed
    # remove x and y labels
    ax.set_xlabel('')
    ax.set_ylabel('')
    y_labels = ax.get_yticklabels()
    ax.set_yticklabels(y_labels, rotation=0)
    # ax.set_title("Pfam domain distribution for top 7 protein families filtered")
    plt.savefig(out_file, format="png", bbox_inches='tight', dpi=800)
    return plt.show()


def dict_count_pfam(df, sp) -> pd.DataFrame:
    # Filter by species and domain type 'Pfam'
    df_ = df.groupby(["sp"]).get_group(sp)
    df_ = df_[df_["db"] == "Pfam"]

    # Count occurrences of PFAM domains
    count = df_.groupby(["family", "db_acc"]).size().reset_index(name=sp)

    # Sort the counts DataFrame first by 'family' in the predefined order, then by count
    sorted_count = count.sort_values(by=['family', sp], ascending=False)

    return sorted_count, df_


def filter_and_sort_families(df, num_entries, family_order, observed_setting=True):
    # Ensure you're working with a copy to avoid setting values on a slice of the original DataFrame
    filtered_df = df.groupby("family", observed=observed_setting).head(num_entries).copy()

    # Convert the 'family' column to a categorical type with the specified order
    filtered_df['family'] = pd.Categorical(filtered_df['family'], categories=family_order, ordered=True)

    # Sort the DataFrame by 'family' according to the predefined order, and then by 'HIN' if needed
    sorted_df = filtered_df.sort_values(by=['family', 'HIN'], ascending=[True, False])
    # sorted_df = filtered_df.sort_values(by=['family', 'trepo'], ascending=[True,False])

    return sorted_df


# Load data
df_ipr = pd.read_csv(ann_ipr_file, header="infer", sep="\t")[["id", "db_acc", "db"]]
df_7sf = pd.read_csv(top7sf_file, header="infer", sep="\t")[["id", "family", "sp"]]
df = pd.merge(df_ipr, df_7sf, how="inner", on="id")
# Initialize count DataFrame with 'HIN' species
count, df_pfam = dict_count_pfam(df, "HIN")
# count, df_pfam = dict_count_pfam(df, "trepo" )

species = ['HIN', 'trepo', 'spiro', 'wb', 'muris', 'carpe', 'kbiala']
for sp in species:
    count = pd.merge(count, dict_count_pfam(df, sp)[0], how="outer")  # get count for pfam domains for each species
    df_pfam = pd.merge(df_pfam, dict_count_pfam(df, sp)[1], how="outer")

counts_top5 = filter_and_sort_families(count, 5, ['LRR', 'HB', 'CP', 'PK', 'CRP', 'WD40', 'ankyrin'])
counts_top10 = filter_and_sort_families(count, 10, ['LRR', 'HB', 'CP', 'PK', 'CRP', 'WD40', 'ankyrin'])
counts_top20 = filter_and_sort_families(count, 20, ['LRR', 'HB', 'CP', 'PK', 'CRP', 'WD40', 'ankyrin'])

plot_heatmap(counts_top10, out_file2)
plot_heatmap(counts_top5, out_file3)
# plot_heatmap( counts_top20, out_file4 )
