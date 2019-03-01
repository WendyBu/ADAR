import pandas as pd
import glob
import os.path
pd.set_option("display.max_column", 100)



def add_identifier(file):
    df = pd.read_csv(file, sep="\t")
    df["ID"] = df["Region"] + "_" + df["Position"].astype(str)
    df.set_index("ID", inplace=True)
    df.to_csv(file, sep="\t")
    pass


def add_ref(file):
    df = pd.read_csv(file, sep="\t")
    ref_df = pd.read_csv("../database/hg19_AG_editing_reference.txt", sep="\t")
    new_df = df.merge(ref_df, left_on="ID", right_on="identifier", how="inner")
    outputF = "../AGTC_ref/" + os.path.basename(file)
    new_df.to_csv(outputF, sep="\t", index=False)
    pass







# step1
# for file in glob.glob("../combineAG&TC/*"):
#     add_identifier(file)
# print "add identifier done!"

#step2
for file in glob.glob("../combineAG&TC/*"):
    add_ref(file)
print "add ref done!"
