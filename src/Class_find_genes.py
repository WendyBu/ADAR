import pandas as pd
from matplotlib import pyplot as plt
import scipy.stats as stats
import numpy as np
pd.set_option("display.max_column", 100)


def find_all_genes(df, refList):
    print len(refList)
    FinalList = []
    for item in refList:
        if item in df.index.tolist():
            FinalList.append(item)
    print len(FinalList)
    sampleTable = df.loc[FinalList, :]
    sampleTable.sort_values("diff", ascending=False, inplace=True)
    sampleTable.fillna("0", inplace=True)
    return sampleTable


def add_annot(df):
    ref = pd.read_csv("../database/hg19_AG_editing_reference.txt", sep="\t", index_col=0)
    ref = ref[["gene", "annot1", "annot2"]]
    return df.join(ref, how="left")


def filter_3UTR(df):

    pass




def filter_5UTR(df):

    pass



def filter_intron(df):

    pass




