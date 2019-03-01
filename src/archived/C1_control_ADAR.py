## Input file:   clean_AGTC
## already referenced, already combined AG&TC,  cleaned up with AG+ and TC-


import pandas as pd
import math
import glob
import os.path
from matplotlib import pyplot as plt
import scipy.stats as stats
import numpy as np
from scipy import ndimage
from matplotlib_venn import venn2
from matplotlib_venn import venn3, venn3_circles
import collections
from scipy.stats import ks_2samp
pd.set_option("display.max_column", 100)


def get_frequency_count_table(inputFile):
    df = pd.read_csv(inputFile, sep="\t", index_col=1)
    df_freq = df[["BaseCount[A,C,G,T]", "Frequency"]]
    return df_freq


def get_frequency_table(inputFile):
    df = pd.read_csv(inputFile, sep="\t", index_col=1)
    df_freq = df[["Frequency"]]
    return df_freq


def combine_freq_tables_cont_ADAR(geneName):
    ADAR = pd.DataFrame()
    A1 = get_frequency_table("../clean_AGTC/A1.txt")
    A2 = get_frequency_table("../clean_AGTC/A2.txt")
    control = A1.join(A2, how="outer", lsuffix="_A1", rsuffix="_A2")
    if geneName == "ADAR1":
        A7 = get_frequency_table("../clean_AGTC/A7.txt")
        A8 = get_frequency_table("../clean_AGTC/A8.txt")
        A9 = get_frequency_table("../clean_AGTC/A9.txt")
        A10 = get_frequency_table("../clean_AGTC/A10.txt")
        A7_A8 = A7.join(A8, how="outer", lsuffix="_A7", rsuffix="_A8")
        A9_A10 = A9.join(A10, how="outer", lsuffix="_A9", rsuffix="_A10")
        ADAR = A7_A8.join(A9_A10, how="outer")
    elif geneName == "ADAR2":
        A11 = get_frequency_table("../clean_AGTC/A11.txt")
        A12 = get_frequency_table("../clean_AGTC/A12.txt")
        ADAR = A11.join(A12, how="outer", lsuffix ="_A11", rsuffix= "_A12")
    else:
        print "No ADAR file was found!"
    control_ADAR = control.join(ADAR, how="outer")
    return control_ADAR


# def clean_list(s):
#     s_list = s.tolist()
#     s_list_clean = [x for x in s_list if str(x)!= "nan"]
#     return s_list_clean
#
#
# def histogram_frequency_distribution(df):
#     ctr1 = clean_list(df.Frequency_A1)
#     ctr2 = clean_list(df.Frequency_A2)
#     ADAR_1 = clean_list(df.Frequency_A7)
#     ADAR_2 = clean_list(df.Frequency_A8)
#     ADAR_3 = clean_list(df.Frequency_A9)
#     ADAR_4 = clean_list(df.Frequency_A10)
#
#     # print ks_2samp(ctr1, ctr2)
#     # print ks_2samp(ADAR_1,ADAR_2)
#     # print ks_2samp(ADAR_3, ADAR_4)  # ***
#
#     ctr = ctr1+ctr2
#     ADAR_sh1 = ADAR_1 + ADAR_2
#     ADAR_sh2 = ADAR_3 + ADAR_4
#     ADAR = ADAR_sh1 + ADAR_sh2
#
#     ctrx, bins = np.histogram(ctr, bins=49)
#
#     density_ctr = stats.gaussian_kde(ctr)
#     density_adar = stats.gaussian_kde(ADAR_1)
#
#     n, x, _ = plt.hist(ctr, bins=50, alpha=0.5, label="Control", density=True,  fill=False, edgecolor="black")
#     plt.plot(x, density_ctr(x))
#     m, y, _ = plt.hist(ADAR_1, bins=50, alpha=0.5, label="ADAR", density=True, fill=False, edgecolor="green")
#     plt.plot(y, density_adar(y))
#     # plt.hist(ADAR_sh2, bins=100, alpha=0.5, label="ADAR_sh2", density=True, fill=False, edgecolor="green")
#     plt.legend(loc='upper right')
#     plt.show()
#
#     diff = n - m
#     print diff, len(diff)
#     print bins, len(bins)
#     plt.plot(bins, diff)
#     plt.show()
#
#     # print ks_2samp(ctr, ADAR_sh1)  # ***
#     # print ks_2samp(ctr, ADAR_sh2)
#     # print ks_2samp(ADAR_sh1, ADAR_sh2) # ***
#     # print ks_2samp(ctr, ADAR)  #***
#     pass


def increased_ADAR_sites(df, cutoff = 0.1):
    """
    CREATE the ARDAR editing pool for this project
    input: clean_AGTC files: A1, A2, A7, A8, A9, A10
    Method: get the sites existed in control but not ADAR1
    get the sites frequency high in control but lower in ADAR2
    :return:
    """
    df["control_freq_mean"] = df[["Frequency_A1", "Frequency_A2"]].mean(axis=1, skipna=True)
    # df["ADAR_freq_mean"] = df[["Frequency_A7", "Frequency_A8", "Frequency_A9", "Frequency_A10"]].mean(axis=1, skipna=True)
    df["ADAR_freq_mean"] = (df.ix[:, 2:]).mean(axis=1, skipna=True)
    df["ADAR_diff_ctr"] = df.ADAR_freq_mean - df.control_freq_mean   # increase editing sites
    diff_df  =  df[df["ADAR_diff_ctr"] >= cutoff]
    all_adar_sites_list = []
    length = 0
    if diff_df.shape[0] > 0:
        # x = diff_df.control_freq_mean
        # y = diff_df.ADAR_freq_mean
        # mystat, pvalue = stats.wilcoxon(x, y,  zero_method='wilcox', correction=False)
        # mystat, pvalue = ks_2samp(x, y)
        all_adar_sites_list = diff_df.index.tolist()
        length = len(all_adar_sites_list)
        # print "cutoff ", cutoff
        # print "num of adar sites all ", length  # 4394
        # print "pvalue", pvalue
    return all_adar_sites_list, length


def define_ADAR_sites(df, cutoff = 0.05):
    """
    CREATE the ARDAR editing pool for this project
    input: clean_AGTC files: A1, A2, A7, A8, A9, A10
    Method: get the sites existed in control but not ADAR1
    get the sites frequency high in control but lower in ADAR2
    :return:
    """
    df["control_freq_mean"] = df[["Frequency_A1", "Frequency_A2"]].mean(axis=1, skipna=True)
    # df["ADAR_freq_mean"] = df[["Frequency_A7", "Frequency_A8", "Frequency_A9", "Frequency_A10"]].mean(axis=1, skipna=True)
    df["ADAR_freq_mean"] = (df.ix[:, 2:]).mean(axis=1,skipna=True)
    df["ADAR_diff_ctr"] = df.control_freq_mean - df.ADAR_freq_mean   # increase editing sites
    diff_df  =  df[df["ADAR_diff_ctr"] >= cutoff]
    all_adar_sites_list = []
    length = 0
    if diff_df.shape[0] > 0:
        # x = diff_df.control_freq_mean
        # y = diff_df.ADAR_freq_mean
        # mystat, pvalue = stats.wilcoxon(x, y,  zero_method='wilcox', correction=False)
        # mystat, pvalue = ks_2samp(x, y)
        all_adar_sites_list = diff_df.index.tolist()
        length = len(all_adar_sites_list)
        # print cutoff, length
    return all_adar_sites_list, length


def add_annotation(df):
    ref = pd.read_csv("../database/hg19_AG_editing_reference.txt", sep="\t", index_col=0)
    ref = ref[["gene", "annot1", "annot2"]]
    newdf = df.join(ref, how="left")
    return newdf


def EZH2_sites_inADAR_pool(all_adar_sites, ezh2cutoff):
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
    EZH2_inADAR_Pool["EZH2diffCont"] = EZH2_inADAR_Pool["EZH2_freq_mean"]-EZH2_inADAR_Pool["control_freq_mean"] # full table for EZH2 in adar pool
    # print "EZH2_inADAR_Pool How many ", EZH2_inADAR_Pool.shape[0]
    EZH2_inADAR_Pool = EZH2_inADAR_Pool[(EZH2_inADAR_Pool.EZH2_freq_mean > 0) & (EZH2_inADAR_Pool.control_freq_mean > 0)]
    # print "non 0 EZH2, control in ADAR pool", EZH2_inADAR_Pool.shape
    EZH2_overEditing = EZH2_inADAR_Pool[EZH2_inADAR_Pool.EZH2diffCont > ezh2cutoff]
    EZH2_underEditing = EZH2_inADAR_Pool[EZH2_inADAR_Pool.EZH2diffCont < (-ezh2cutoff)]
    # print "Ezh2 underEditing", EZH2_underEditing.shape   #707
    # print "ezh2 cross adar table ", EZH2_inADAR_Pool.shape  # 3769
    # print "ezh2 overediting ", EZH2_overEditing.shape  #251
    # print "unaffected", EZH2_inADAR_Pool.shape[0] - EZH2_overEditing.shape[0] - EZH2_underEditing.shape[0]  #2811

    EZH2_inADAR_Pool = EZH2_inADAR_Pool.round(3)
    EZH2_underEditing = add_annotation(EZH2_underEditing)
    EZH2_overEditing = add_annotation(EZH2_overEditing)
    EZH2_inADAR_Pool = add_annotation(EZH2_inADAR_Pool)


    # EZH2_overEditing.to_csv("../editingPool2/overEditing.xls", sep="\t")
    # EZH2_underEditing.to_csv("../editingPool2/underEditing.xls", sep="\t")
    # EZH2_inADAR_Pool.to_csv("../editingPool2/EZH2inADARPool.xls", sep="\t")
    return EZH2_inADAR_Pool, EZH2_overEditing, EZH2_underEditing


def getcutoff_ADAR(df):
    """
    decreased adar sites define
    :param df:
    :return:
    """
    n = np.linspace(0,1.0, num=101, endpoint=True)
    cutoff_numOfSites = {}
    for cutoff in n:
        all_adar_sites, numberOfSites= define_ADAR_sites(df, cutoff)  # a list
        cutoff_num_pair = {cutoff: numberOfSites}
        cutoff_numOfSites.update(cutoff_num_pair)
    new_df = pd.DataFrame.from_dict(cutoff_numOfSites, orient='index')
    return new_df


def increased_ADAR(df):
    """
    the sites increased editing freq in ADAR shRNA samples
    :param df:
    :return:
    """
    n = np.linspace(0, 1.0, num=101, endpoint=True)
    cutoff_numOfSites = {}
    for cutoff in n:
        all_adar_sites, numberOfSites = increased_ADAR_sites(df, cutoff)  # a list
        cutoff_num_pair = {cutoff: numberOfSites}
        cutoff_numOfSites.update(cutoff_num_pair)
    new_df = pd.DataFrame.from_dict(cutoff_numOfSites, orient='index')
    return new_df



def plot_wired(ctr_ADAR_freq_df):
    decreased_table = getcutoff_ADAR(ctr_ADAR_freq_df)
    increased_table = increased_ADAR(ctr_ADAR_freq_df)
    decreased_number = decreased_table.shape[0]
    increased_number = increased_table.shape[0]
    print decreased_table.head()
    df = decreased_table.join(increased_table, lsuffix="_decreased", rsuffix="_increase")
    df.sort_index(ascending=True, inplace=True)
    s = np.linspace(0, 100, num=101, endpoint=True)
    s = np.asanyarray(s)
    df["order"] = s
    df["decrease_percentage"] =  df["0_decreased"]/decreased_number
    df["increase_percentage"] = df["0_increase"]/increased_number
    # df.to_csv("../editingPool2/table_increase_decrease_adar.xls", sep="\t")
    x = df["increase_percentage"]
    y = df["decrease_percentage"]

    plt.plot(x, y)
    plt.show()
    # fig, ax1 = plt.subplots()
    # ax1.plot(x, y)
    # plt.gca().invert_xaxis()
    # plt.gca().invert_yaxis()
    # ax2 = ax1.twinx()
    # ax3 = ax1.twiny()
    # cutoff_scale = df.loc[:, "order"].tolist()
    # cutoff_scale = cutoff_scale[::-1]
    # ax2.tick_params(cutoff_scale)
    # ax3.tick_params(cutoff_scale)
    # plt.show()
    return df


def find_shortest_distance(df):
    df["y"] = (df["decrease_percentage"]-1)**2.0
    df["x"] = df["increase_percentage"] ** 2.0
    df["dist"] = np.sqrt(df["x"] + df["y"])
    shortest_dist = df.dist.min()
    optimize_cutoff = df[["dist"]].idxmin()
    # df.to_csv("../editingPool2/distance.xls", sep="\t")
    cutoffs = df.index.tolist()
    distance = df.dist.tolist()
    plt.plot(cutoffs, distance)
    plt.show()
    return optimize_cutoff


def main():
    ################### Get the frequency Table###############################
    ctr_ADAR_freq_df = combine_freq_tables_cont_ADAR("ADAR2")  # or ADAR1

    ####################Look for ADAR cutoff###################################
    ROC_data = plot_wired(ctr_ADAR_freq_df)
    ADAR_cutoffs = find_shortest_distance(ROC_data)
    print ADAR_cutoffs










    #histogram_frequency_distribution(ctr_ADAR_freq_df)   # the frequency distribution map of three samples

    # plot the specfic figures to look for cutoff.


    ############## all adar changed sites##########################################
    # all_decreased_size, totalNum_decrease = define_ADAR_sites(ctr_ADAR_freq_df, cutoff=0.00)
    # all_increased_size, totalNum_increase = increased_ADAR_sites(ctr_ADAR_freq_df, cutoff=0.00)
    # increased_sites = ctr_ADAR_freq_df.loc[all_increased_size, :]
    # decreased_sites = ctr_ADAR_freq_df.loc[all_decreased_size, :]
    # decreased_sites.to_csv("../editingPool2/all_decreased_sites.xls", sep="\t")
    # increased_sites.to_csv("../editingPool2/all_increased_sites.xls", sep="\t")
    # print "total number of decrease in adar sh ", totalNum_decrease
    # print "total number of increase in adar sh ", totalNum_increase  ##8983, 5394.


    # step1 : loop cutoff and decide what cutoff to use to define ADAR sites
    # cutoff_numberSites = getcutoff_ADAR()

    # step2 : use the decided cutoff to generate ADAR sites
    # all_adar_sites, num_sites = define_ADAR_sites(ctr_ADAR_freq_df, 0.09)
    # print num_sites
    # all_adar_sites = list(set(all_adar_sites))
    #
    # # # step 3 : cross with EZH2 sites
    # ezh2_df, up_df, down_df = EZH2_sites_inADAR_pool(all_adar_sites, 0.10)
    # # print ezh2_df.groupby(["annot1"]).size()
    # print up_df.groupby(["annot1"]).size()
    # print down_df.groupby(["annot1"]).size()
    pass


if __name__ == "__main__":
    main()



# work for both ADAR1 and ADAR2
# input has to be clear ADAR1 or ADAR2