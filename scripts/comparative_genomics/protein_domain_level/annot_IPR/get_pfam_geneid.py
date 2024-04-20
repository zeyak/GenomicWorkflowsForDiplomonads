import pandas as pd

"""
This script;
gets the gene_ids with the given PFAM domains
"""
pfam= "PF12799"

top7sf = "output/interproscan/processed_data/ann_ipr_cat_top7sf.csv" #get pfam domains
HIN_interpro = "output/interproscan/processed_data/HIN.csv" #get pfam domains

out_file = f"""output/interproscan/processed_data/{pfam}_geneid.csv""" # get gene_ids with pfam domains


def get_pfam_geneid(df, out_file, pfam):
    df = df[df["db_acc"] == pfam]
    df = df[["id", "db_acc"]].drop_duplicates(subset="id").reset_index(drop=True)
    print(pfam, len(df))
    df["id"].to_csv(out_file, index=False, sep="\t", header=None)
    return df

df_top7 = pd.read_csv(top7sf, header="infer", sep="\t")["id"] #get gene ids from top 7 superfamily
df_HIN = pd.read_csv(HIN_interpro, header="infer", sep="\t") #ger gene ids from HIN

pfam_HIN = pd.merge(df_top7, df_HIN, how="inner")[["id", "db_acc"]].drop_duplicates() #get top7 superfamily gene ids with pfam domains

df_PF12799 = get_pfam_geneid(pfam_HIN, out_file,pfam ) #write gene ids with PF12799
#df_PF13855 = get_pfam_geneid(HIN_interpro, pfam_geneid_HIN, "PF13855")
#df_PF13516 = get_pfam_geneid(HIN_interpro, pfam_geneid_HIN, "PF13516")