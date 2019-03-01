## Input file:   clean_AGTC
## already referenced, already combined AG&TC,  cleaned up with AG+ and TC-


import pandas as pd
import glob
import os.path
from matplotlib import pyplot as plt
import numpy as np
from matplotlib_venn import venn2
from matplotlib_venn import venn3, venn3_circles
import collections
pd.set_option("display.max_column", 100)


def get_frequency_table(inputFile):
    df = pd.read_csv(inputFile, sep="\t", index_col=1)
    df_freq = df[["BaseCount[A,C,G,T]", "Frequency", "gene"]]
    return df_freq


def ADAR_only_sites(cutoff=0.67):
    """
    CREATE the ARDAR editing pool for this project
    input: clean_AGTC files: A1, A2, A7, A8, A9, A10
    Method: get the sites existed in control but not ADAR1
    get the sites frequency high in control but lower in ADAR2
    :return:
    """

    A1 = get_frequency_table("../clean_AGTC/A1.txt")
    A2 = get_frequency_table("../clean_AGTC/A2.txt")
    A7 = get_frequency_table("../clean_AGTC/A7.txt")
    A8 = get_frequency_table("../clean_AGTC/A8.txt")
    A9 = get_frequency_table("../clean_AGTC/A9.txt")
    A10 = get_frequency_table("../clean_AGTC/A10.txt")
    control = A1.join(A2, how = "outer", lsuffix="_A1", rsuffix="_A2")
    A7_A8 = A7.join(A8, how="outer", lsuffix = "_A7", rsuffix = "_A8")
    A9_A10 = A9.join(A10, how="outer", lsuffix="_A9", rsuffix="_A10")
    ADAR = A7_A8.join(A9_A10, how="outer")
    control_ADAR = control.join(ADAR, how="outer")
    control_ADAR["control_freq_mean"] = control_ADAR[["Frequency_A1", "Frequency_A2"]].mean(axis=1, skipna=True)
    control_ADAR["ADAR_freq_mean"] = control_ADAR[["Frequency_A7", "Frequency_A8", "Frequency_A9", "Frequency_A10"]].mean(axis=1, skipna=True)
    control_ADAR["ADARvsCont"] = (control_ADAR.ADAR_freq_mean+0.01) / (control_ADAR.control_freq_mean+0.01)
    ADAR_diff_cont = control_ADAR[control_ADAR.ADARvsCont < cutoff]
    all_adar_sites_list = ADAR_diff_cont.index.tolist()
    length = len(all_adar_sites_list)
    print "cutoff ", cutoff
    print "num of adar sites all ", length  # 4394
    return all_adar_sites_list, length


def add_annotation(df):
    ref = pd.read_csv("../database/hg19_AG_editing_reference.txt", sep="\t", index_col=0)
    ref = ref[["gene", "annot1", "annot2"]]
    newdf = df.join(ref, how="left")
    return newdf


def EZH2_sites_inADAR_pool(all_adar_sites):
    """
    adar sites pool
    :return: under_editing; over_editing; No_change_in_EZH2vsControl
    """
    A1 = get_frequency_table("../clean_AGTC/A1.txt")
    A2 = get_frequency_table("../clean_AGTC/A2.txt")
    A3 = get_frequency_table("../clean_AGTC/A3.txt")
    A4 = get_frequency_table("../clean_AGTC/A4.txt")
    A5 = get_frequency_table("../clean_AGTC/A5.txt")
    A6 = get_frequency_table("../clean_AGTC/A6.txt")
    control = A1.join(A2, how="outer", lsuffix="_A1", rsuffix="_A2")
    A3_A4 = A3.join(A4, how="outer", lsuffix="_A3", rsuffix="_A4")
    A5_A6 = A5.join(A6, how="outer", lsuffix="_A5", rsuffix="_A6")
    EZH2 = A3_A4.join(A5_A6, how="outer")
    control_EZH2 = control.join(EZH2, how="outer")
    control_EZH2["control_freq_mean"] = control_EZH2[["Frequency_A1", "Frequency_A2"]].mean(axis=1, skipna=True)
    control_EZH2["EZH2_freq_mean"] = control_EZH2[["Frequency_A3", "Frequency_A4", "Frequency_A5", "Frequency_A6"]].mean(axis=1, skipna=True)

    EZH2_inADAR_Pool = control_EZH2.loc[all_adar_sites, :]   # all EZH2 and control editing sites
    EZH2_inADAR_Pool["EZH2vsContRatio"] = (EZH2_inADAR_Pool["EZH2_freq_mean"]+0.01)/(EZH2_inADAR_Pool["control_freq_mean"]+0.01)  # full table for EZH2 in adar pool
    print "EZH2_inADAR_Pool How many ", EZH2_inADAR_Pool.shape[0]
    EZH2_inADAR_Pool = EZH2_inADAR_Pool[(EZH2_inADAR_Pool.EZH2_freq_mean > 0) & (EZH2_inADAR_Pool.control_freq_mean > 0)]
    print "non 0 EZH2, control in ADAR pool", EZH2_inADAR_Pool.shape
    EZH2_overEditing = EZH2_inADAR_Pool[EZH2_inADAR_Pool.EZH2vsContRatio > 1.5]
    EZH2_underEditing = EZH2_inADAR_Pool[EZH2_inADAR_Pool.EZH2vsContRatio < 0.67]
    print "Ezh2 underEditing", EZH2_underEditing.shape   #707
    print "ezh2 cross adar table ", EZH2_inADAR_Pool.shape  # 3769
    print "ezh2 overediting ", EZH2_overEditing.shape  #251
    print "unaffected", EZH2_inADAR_Pool.shape[0] - EZH2_overEditing.shape[0] - EZH2_underEditing.shape[0]  #2811


    EZH2_underEditing = add_annotation(EZH2_underEditing)
    EZH2_overEditing = add_annotation(EZH2_overEditing)
    EZH2_inADAR_Pool = add_annotation(EZH2_inADAR_Pool)

    EZH2_overEditing.to_csv("../editingPool/overEditing.xls", sep="\t")
    EZH2_underEditing.to_csv("../editingPool/underEditing.xls", sep="\t")
    EZH2_inADAR_Pool.to_csv("../editingPool/EZH2inADARPool.xls", sep="\t")
    pass


def getcutoff_ADAR():
    getcutoff_ADAR()
    n = [1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
    cutoff_numOfSites = {}
    for cutoff in n:
        all_adar_sites, numberOfSites = ADAR_only_sites(cutoff)  # a list
        cutoff_num_pair = {cutoff: numberOfSites}
        cutoff_numOfSites.update(cutoff_num_pair)
    cutoff_numOfSites = collections.OrderedDict(sorted(cutoff_numOfSites.items()))
    return cutoff_numOfSites


def main():
    # step1 : loop cutoff and decide what cutoff to use to define ADAR sites
    # cutoff_numberSites = getcutoff_ADAR()
    # step2 : use the decided cutoff to generate ADAR sites
    all_adar_sites = ADAR_only_sites(cutoff = 0.67)
    # step 3 : cross with EZH2 sites
    EZH2_sites_inADAR_pool(all_adar_sites)
    pass


if __name__ == "__main__":
    main()



