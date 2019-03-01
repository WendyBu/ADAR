"""
input: organized df
eg
                   A1  A2   A7   A8    A9   A10  standardMean  expMean   diff
ID
chr10_101144703   NaN NaN  0.4  NaN   NaN   NaN           NaN    0.400    NaN
chr10_101159356  0.12 NaN  NaN  NaN  0.48  0.17          0.12    0.325  0.205
chr10_101159524  0.11 NaN  NaN  NaN  0.10   NaN          0.11    0.100 -0.010
chr10_101159528   NaN NaN  NaN  NaN  0.10   NaN           NaN    0.100    NaN
chr10_101160526   NaN NaN  NaN  0.2   NaN   NaN           NaN    0.200    NaN

output: ROC table
"""

import pandas as pd
from matplotlib import pyplot as plt
import scipy.stats as stats
import numpy as np
pd.set_option("display.max_column", 100)


def getDiff_df(df, change=""):
    """
    :param df: input df: with full frequency; use by diffAll()
    :param change: "increase" or "decrease"
    :return: increased df or decreased df
    """
    df_diff, change_num = pd.DataFrame(), 0
    if change == "increase":
        df_diff = df[df["diff"] > 0]
        change_num = df_diff.shape[0]
        print "The increased site number is ", change_num
    elif change == "decrease":
        df_diff = df[df["diff"] < 0]
        change_num = df_diff.shape[0]
        print "The decreased site number is ", change_num
    else:
        print "cannot find change information!"
    return df_diff, change_num


def diffAll(df):
    """
    :param df: frequency table; Work with getDiff_df()
    :return: table with decrease; table with increase; and numbers of sites
    """
    decreased_table, decreased_number = getDiff_df(df, change="decrease")
    increased_table, increased_number = getDiff_df(df, change="increase")
    decrease_df = (-decreased_table["diff"]).sort_values(ascending = False)  #change all the negative value to positive
    increase_df = increased_table["diff"].sort_values(ascending = False)
    return decrease_df, increase_df, decreased_number, increased_number


def define_sites(df, cutoff):
    diff_df = df[df.values > cutoff]  # the value is bigger than cutoff(only one col in the input df
    all_adar_sites_list = diff_df.index.tolist()
    length = len(all_adar_sites_list)
    return all_adar_sites_list, length


def loopCutoff(df):
    """
    using cutoff 0-1.0(100 points), accumulate the sites number for each cutoff
    :param df: decreased df or increased df, work with define_sites()
    :return: df, index is cutoff, value is the accumulates number of sites
    """
    n = np.linspace(0,1.0, num=101, endpoint=True)  #loop 100 times
    cutoff_numOfSites = {}
    for cutoff in n:
        all_sites, numberOfSites= define_sites(df, cutoff)  # a list
        cutoff_num_pair = {cutoff: numberOfSites}
        cutoff_numOfSites.update(cutoff_num_pair)
    newdf = pd.DataFrame(list(cutoff_numOfSites.items()), columns=["cutoff", "SitesNumber"])
    newdf.sort_values("cutoff", ascending=True, inplace=True)
    newdf.set_index("cutoff", inplace=True)
    return newdf


def plotROC(df):
    y = df.decrease_percentage
    x = df.increase_percentage
    plt.plot(x, y)
    plt.show()
    pass


def find_shortest_distance(df):
    df["y"] = (df["decrease_percentage"]-1)**2.0
    df["x"] = df["increase_percentage"] ** 2.0
    df["dist"] = np.sqrt(df["x"] + df["y"])
    optimize_cutoff = df[["dist"]].idxmin()
    cutoffs = df.index.tolist()
    distance = df.dist.tolist()
    plt.plot(cutoffs, distance)
    plt.show()
    return optimize_cutoff.values[0]


def ROC(df):
    decreased_table, increased_table, decreased_number, increased_number = diffAll(df)
    decreaseAccuSitesCutoff = loopCutoff(decreased_table)
    increaseAccuSitesCutoff = loopCutoff(increased_table)
    IncreaseDecreaseNumSites = decreaseAccuSitesCutoff.join(increaseAccuSitesCutoff, lsuffix="_d", rsuffix="_i")
    IncreaseDecreaseNumSites["decrease_percentage"] = IncreaseDecreaseNumSites["SitesNumber_d"]/decreased_number
    IncreaseDecreaseNumSites["increase_percentage"] = IncreaseDecreaseNumSites["SitesNumber_i"]/increased_number
    plotROC(IncreaseDecreaseNumSites)
    final_cutoff = find_shortest_distance(IncreaseDecreaseNumSites)
    siteList, length = define_sites(decreased_table, final_cutoff)  # only care about sites decreased compared to control
    return final_cutoff, siteList, length


"""
table: IncreaseDecreaseNumSites
        SitesNumber_d  SitesNumber_i  decrease_percentage  increase_percentage
cutoff                                                                        
0.00             8574           4985             1.000000             1.000000
0.01             8252           4699             0.962445             0.942628
0.02             7781           4243             0.907511             0.851153
0.03             7286           3804             0.849778             0.763089
0.04             6926           3525             0.807791             0.707121

"""




