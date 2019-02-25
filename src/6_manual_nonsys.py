import pandas as pd

df = pd.read_table("../NonSys/new/control.txt", sep="\t")
print df
list1 = df.iloc[:, 0].tolist()
list2 = df.iloc[:, 1].tolist()
common_gene_control = set(list1).intersection(set(list2))

