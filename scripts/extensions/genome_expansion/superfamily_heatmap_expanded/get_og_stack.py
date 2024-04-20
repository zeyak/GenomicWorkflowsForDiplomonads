import pandas as pd

try:
    og_file = snakemake.input.og
    out = snakemake.output[0]
except NameError:
    # testing
    og_file = "output/1_orthofinder/Results_Oct17_2/Orthogroups/Orthogroups.txt"
    out_file = "data/orthogroups/og_stack.csv"

df = pd.read_csv(og_file, header=None, dtype=str, delim_whitespace=True)
df = df.apply(lambda x: x.str.replace(":", ""))
df = df.rename(columns={0: "OG"}).set_index("OG")
df = df.stack().reset_index().drop(columns=["level_1"])
df = df.rename(columns={0: "id"})


def convert_name(name: str) -> str:
    """
    Convert names.
    @param df:
    @return:
    """
    pairs = {
        'HIN': 'HIN',
        'TPC1': 'trepo',
        'SS': 'spiro',
        'GL': 'wb',
        'GMRT': 'muris',
        'KAG': 'carpe',
        'GIQ': 'kbiala',
        'GCA': 'kbiala',

    }

    for key, value in pairs.items():
        if name.startswith(key):
            return value


df['sp'] = df.id.apply(convert_name)

df.to_csv(out_file, sep="\t", index=False)
