import pandas as pd
import glob
import os.path
pd.set_option("display.max_column", 100)

### add the number of total copy
for file in glob.glob("../NonSys/combine/*"):
    df = pd.read_csv(file, sep="\t")
    df["ACGT"] = df.loc[:, "BaseCount[A,C,G,T]"].str.extract(r'(\d+\S\s\d+\S\s\d+\S\s\d+)')
    df[["A", "C", "G","T"]] = df.loc[:, "ACGT"].str.split(",", expand=True)
    df.A = pd.to_numeric(df.A)
    df.C = pd.to_numeric(df.C)
    df.G = pd.to_numeric(df.G)
    df["T"] = pd.to_numeric(df["T"])
    df["total_copy"] = df.loc[:, ["A", "C", "G", "T"]].sum(axis=1)
    del df["ACGT"]
    df.to_csv(file + ".xls", sep="\t")

