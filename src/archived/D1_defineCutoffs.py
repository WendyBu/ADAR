## Input file:   clean_AGTC
## already referenced, already combined AG&TC,  cleaned up with AG+ and TC-


import pandas as pd
from matplotlib import pyplot as plt
import scipy.stats as stats
import numpy as np
pd.set_option("display.max_column", 100)



def add_annotation(df):
    ref = pd.read_csv("../database/hg19_AG_editing_reference.txt", sep="\t", index_col=0)
    ref = ref[["gene", "annot1", "annot2"]]
    newdf = df.join(ref, how="left")
    return newdf








def plot_ROC(df):
    decreased_table = ADAR_loop(df, change="decrease")
    increased_table = ADAR_loop(df, change="increase")
    df_diff = decreased_table.join(increased_table, lsuffix="_decreased", rsuffix="_increase")
    df_diff.sort_index(ascending=True, inplace=True)
    increased_number, decreased_number = increase_decrease_number(df)

    df_diff["decrease_percentage"] =  (df_diff["0_decreased"]/decreased_number)
    df_diff["increase_percentage"] = (df_diff["0_increase"]/increased_number)
    x = df_diff["increase_percentage"]
    y = df_diff["decrease_percentage"]
    plt.plot(x, y)
    plt.show()
    return df_diff


def find_shortest_distance(df):
    df["y"] = (df["decrease_percentage"]-1)**2.0
    df["x"] = df["increase_percentage"] ** 2.0
    df["dist"] = np.sqrt(df["x"] + df["y"])
    optimize_cutoff = df[["dist"]].idxmin()
    cutoffs = df.index.tolist()
    distance = df.dist.tolist()
    plt.plot(cutoffs, distance)
    plt.show()
    return optimize_cutoff


def EZH2_sites_inADAR_pool(all_adar_sites, ezh2cutoff):
    """
    adar sites pool
    :return: under_editing; over_editing; No_change_in_EZH2vsControl
    """
    control = A1.join(A2, how="outer", lsuffix="_A1", rsuffix="_A2")
    A3_A4 = A3.join(A4, how="outer", lsuffix="_A3", rsuffix="_A4")
    A5_A6 = A5.join(A6, how="outer", lsuffix="_A5", rsuffix="_A6")
    EZH2 = A3_A4.join(A5_A6, how="outer")
    control_EZH2 = control.join(EZH2, how="outer")
    control_EZH2["control_freq_mean"] = control_EZH2[["Frequency_A1", "Frequency_A2"]].mean(axis=1, skipna=True)
    control_EZH2["EZH2_freq_mean"] = control_EZH2[["Frequency_A3", "Frequency_A4", "Frequency_A5", "Frequency_A6"]].mean(axis=1, skipna=True)

    EZH2_inADAR_Pool = control_EZH2.loc[all_adar_sites, :]
    EZH2_inADAR_Pool["EZH2diffCont"] = EZH2_inADAR_Pool["EZH2_freq_mean"]-EZH2_inADAR_Pool["control_freq_mean"] # full table for EZH2 in adar pool
    EZH2_inADAR_Pool = EZH2_inADAR_Pool[(EZH2_inADAR_Pool.EZH2_freq_mean > 0) & (EZH2_inADAR_Pool.control_freq_mean > 0)]
    EZH2_overEditing = EZH2_inADAR_Pool[EZH2_inADAR_Pool.EZH2diffCont > ezh2cutoff]
    EZH2_underEditing = EZH2_inADAR_Pool[EZH2_inADAR_Pool.EZH2diffCont < (-ezh2cutoff)]
    EZH2_underEditing = add_annotation(EZH2_underEditing)
    EZH2_overEditing = add_annotation(EZH2_overEditing)
    EZH2_inADAR_Pool = add_annotation(EZH2_inADAR_Pool)
    EZH2_inADAR_Pool.to_csv("../editingPool2/ADAR2_EZH2.xls", sep="\t")
    under = EZH2_underEditing.shape[0]
    over = EZH2_overEditing.shape[0]
    all = EZH2_inADAR_Pool.shape[0]
    unaffected = all - under - over
    print "EZH2 cutoff", ezh2cutoff
    print "EZH2 under Editing", under
    print "EZH2 over Editing", over
    print "EZH2 unaffected", unaffected
    # EZH2_overEditing.to_csv("../editingPool2/overEditing.xls", sep="\t")
    # EZH2_underEditing.to_csv("../editingPool2/underEditing.xls", sep="\t")
    # EZH2_inADAR_Pool.to_csv("../editingPool2/EZH2inADARPool.xls", sep="\t")
    return  under, over, unaffected


def main():

    ROC_data = plot_ROC(ctr_ADAR_freq_df)
    ADAR_cutoffs = find_shortest_distance(ROC_data)
    ADAR_cutoff_value = ADAR_cutoffs.values[0]
    print "ADAR cuttoff", ADAR_cutoff_value
    ######################EZH2_sites_inADAR_pool###############
    all_adar_sites, num_adar_sites = define_ADAR_sites(ctr_ADAR_freq_df, cutoff = ADAR_cutoff_value, change="decrease")
    # ezh2cutoffs = [0.05, 0.10, 0.15, 0.20]
    # ezh2_cutoff_numbers = pd.DataFrame()
    # for ezh2cutoff in ezh2cutoffs:
    #     under, over, unaffected = EZH2_sites_inADAR_pool(all_adar_sites, ezh2cutoff)
    #     ezh2_cutoff_numbers = ezh2_cutoff_numbers.append({ezh2cutoff: [under, over, unaffected]})
    # print ezh2_cutoff_numbers
    under, over, unaffected = EZH2_sites_inADAR_pool(all_adar_sites, 0.05)
    pass


if __name__ == "__main__":
    main()



# work for both ADAR1 and ADAR2
# input has to be clear ADAR1 or ADAR2