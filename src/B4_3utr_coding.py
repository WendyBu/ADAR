import pandas as pd
import glob
import os.path
from matplotlib import pyplot as plt
import numpy as np
from matplotlib_venn import venn2
from matplotlib_venn import venn3, venn3_circles
pd.set_option("display.max_column", 100)


def getsubset(file, key):
    df = pd.read_csv(file, sep="\t", index_col=0)
    df_key = df[df.annot1 == key]
    return df_key


def main():
    underEditing_3UTR = getsubset("../editingPool/underEditing.xls", "3UTR")  #70 items
    overEditing_3UTR = getsubset("../editingPool/overEditing.xls", "3UTR")
    underEditing_3UTR.to_csv("../editingPool/underEditing_3UTR.xls",sep="\t")
    overEditing_3UTR.to_csv("../editingPool/overEditing_3UTR.xls",sep="\t")
    all_coding = getsubset("../editingPool/EZH2inADARPool.xls", "Nonsyn")
    all_coding.to_csv("../editingPool/allcoding.xls", sep="\t")
    pass


if __name__ == "__main__":
    main()