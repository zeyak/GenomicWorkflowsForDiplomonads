import pandas as pd

try:
    og_file = snakemake.input.og
    out = snakemake.output[0]
except NameError:
    # testing
    og_ann_lgt = "data/orthogroups/og_ann_lgt.csv"
    og_count_s = "data/orthogroups/og_count_s.csv"

    out_file_csv = "data/og_comb/{}.csv"
    out_file_xlsx= "data/og_comb/{}.xlsx"


def get_og_genes(og_ann_lgt, og_category: pd.DataFrame) -> pd.DataFrame:
    """
    @param og_category: Orthgroups with a species combination.
    @return: Orthogroup with its genes for given category.
    """
    og_ann_lgt_ = pd.read_csv(og_ann_lgt, sep="\t", header="infer")
    df = pd.merge(og_ann_lgt_, og_category)

    return df


df = pd.read_csv(og_count_s, sep="\t", header='infer')


""" OG Combinations """
og_comb = {
    "og_diplo": df[(df.loc[:, ["HIN", "trepo", "spiro", "wb", "muris"]] >= 1).all(1) &
                   (df.loc[:, ["carpe", "kbiala"]] == 0).all(1)],

    "og_hin_trepo" : df[(df.loc[:, ["HIN", "trepo"]] >= 1).all(1) &
                        (df.loc[:, ["spiro", "wb", "muris", "carpe", "kbiala"]] == 0).all(1)],

    "og_hin_trepo" : df[(df.loc[:, ["HIN", "trepo"]] >= 1).all(1) &
                        (df.loc[:, ["spiro", "wb", "muris", "carpe", "kbiala"]] == 0).all(1)],

    "og_fl" : df[(df.loc[:, ["carpe", "kbiala", "HIN", "trepo"]] >= 1).all(1) &
                 (df.loc[:, ["spiro", "wb", "muris"]] == 0).all(1)],

    "og_hin_trepo_kbiala" : df[(df.loc[:, ["HIN", "trepo", "kbiala"]] >= 1).all(1) &
                               (df.loc[:, ["spiro", "wb", "muris", "carpe"]] == 0).all(1)],

    "og_hin_trepo_carpe" : df[(df.loc[:, ["HIN", "trepo", "carpe"]] >= 1).all(1) &
                              (df.loc[:, ["spiro", "wb", "muris", "kbiala"]] == 0).all(1)]
}


dic_groups = {}
for keys, values in og_comb.items():

    dic_groups[keys] = get_og_genes(og_ann_lgt, values).sort_values(by=["LGT", "OG", "Total"])

    dic_groups[keys].to_csv(out_file_csv.format(keys), sep="\t", index=False)
    dic_groups[keys].to_excel(out_file_xlsx.format(keys), index=False)


dic_stats={}
""" og_comb statistics """
for keys, values in og_comb.items():
    dic_stats[keys] = values.agg("sum", axis=0)

