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


