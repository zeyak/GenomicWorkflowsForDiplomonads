import pandas as pd
import glob

"""
This script;
Filters each interproscan output from diff species
Concatanates all the outputs
"""

try:
    og_file = snakemake.input.og
    out = snakemake.output[0]
except NameError:
    # testing
    path = "output/interproscan/raw_data/*.tsv"
    out_file = "output/interproscan/processed_data/ann_ipr_cat.csv"
    list_files = glob.glob(path)


#get species name from the filenames
sp_dic = {}
for element in list_files:
    i = element.split(".")[0]
    i = i.split("/")[-1]
    sp_dic[i] = element


dic_ann= {}
for key, value in sp_dic.items():
    dic_ann[key] = pd.read_csv(value, sep="\t", header=None, names=list(range(0, 15)),
                               engine='python', quoting=3)[[0,3,4,5,11,12]]
    dic_ann[key] = dic_ann[key].dropna().drop_duplicates().rename(
        columns={0: "id", 3:"db", 4:"db_acc", 5:"ann_db", 11:"ipr", 12:"ann_inter"})
    dic_ann[key]["sp"]= key
    dic_ann[key].to_csv(f"output/interproscan/processed_data/{key}.csv", sep="\t", index=False)
    

pd.concat(dic_ann, axis=0).dropna().drop_duplicates().to_csv(out_file, sep="\t", index=False)
