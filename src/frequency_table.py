import Class_combine_freq_table as cmb
import Class_define_cutoff as cf
import Class_find_genes as fg
import pandas as pd


def main():
    df_c_e= cmb.organize_group(standard="control", exp="EZH2")
    df_ADAR = cmb.organize_group("ADAR1", "ADAR2")
    df = df_c_e.join(df_ADAR, lsuffix = "_c_e", rsuffix = "_ADAR", how = "outer")
    df_sample = df[["A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10", "A11", "A12"]]
    df_stat = pd.DataFrame()
    df_stat[["Control", "EZH2", "ADAR1", "ADAR2"]] = df[["standardMean_c_e", "expMean_c_e", "standardMean_ADAR", "expMean_ADAR"]]
    df_sample.fillna("0", inplace=True)
    df_stat.fillna("0", inplace=True)
    df_sample.to_csv("../tempResults/sample_frequency_total_table.xls", sep="\t")
    df_stat.to_csv("../tempResults/sample_frequency_average.xls", sep="\t")
    return df_sample, df_stat


if __name__ == "__main__":
    main()