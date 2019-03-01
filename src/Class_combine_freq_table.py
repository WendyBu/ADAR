"""input: clean_AGTC
12 txt files
return: frequency table
"""

import pandas as pd
import os.path
pd.set_option("display.max_column", 100)


def get_frequency_table(inputFile, file):
    df = pd.read_csv(inputFile, sep="\t", index_col=1)
    df_freq = df[["Frequency"]]
    df_freq.rename(index = str, columns={"Frequency":file}, inplace=True)
    return df_freq


def combine_InGroup(filenames):
    group_df = pd.DataFrame()
    inputDir = "../clean_AGTC"  # the txt file storage dir
    for file in filenames:
        inputFile = os.path.join(inputDir, file + ".txt")
        df = get_frequency_table(inputFile, file)
        if group_df.empty:
            group_df = df
        else:
            group_df = group_df.join(df, how="outer")
    return group_df


def combineGroups(standard_df, exp_df):
    return standard_df.join(exp_df, how="outer")


def computeStats(df, standardSample, expSample):
    df["standardMean"] = df[standardSample].mean(axis=1, skipna=True)
    df["expMean"] = df[expSample].mean(axis=1, skipna=True)
    df["diff"] = df.expMean - df.standardMean
    return df


def organize_group(standard="control", exp = "ADAR1"):
    sample_info = {"control":["A1", "A2"], "EZH2":["A3", "A4", "A5", "A6"], "ADAR1":["A7", "A8", "A9", "A10"], "ADAR2":["A11", "A12"]}
    standardSample = sample_info.get(standard)
    expSample  = sample_info.get(exp)
    standard_df = combine_InGroup(standardSample)   #standard_frequency dataframe
    exp_df = combine_InGroup(expSample)             #exp_frequency dataframe
    standard_exp_df = combineGroups(standard_df, exp_df)
    stat_df = computeStats(standard_exp_df, standardSample, expSample)
    return stat_df



"""
give standard="control", exp = "ADAR1"
return organized df
                  A1  A2   A7   A8    A9   A10  standardMean  expMean   diff
ID                                                                           
chr10_101144703   NaN NaN  0.4  NaN   NaN   NaN           NaN    0.400    NaN
chr10_101159356  0.12 NaN  NaN  NaN  0.48  0.17          0.12    0.325  0.205
chr10_101159524  0.11 NaN  NaN  NaN  0.10   NaN          0.11    0.100 -0.010
chr10_101159528   NaN NaN  NaN  NaN  0.10   NaN           NaN    0.100    NaN
chr10_101160526   NaN NaN  NaN  0.2   NaN   NaN           NaN    0.200    NaN

"""