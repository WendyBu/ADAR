import Class_combine_freq_table as cmb
import Class_define_cutoff as cf
import Class_find_genes as fg
import numpy as np
import scipy.stats as stat
from matplotlib import pyplot as plt


class Three_sample_corr:
    def __init__(self, sample1, sample2, sample3):
        self.sample1 = sample1
        self.sample2 = sample2
        self.sample3 = sample3


    def combine_tables(self):
        df1 = cmb.organize_group(standard="control", exp=self.sample1)
        df2 = cmb.organize_group(standard="control", exp=self.sample2)
        df3 = cmb.organize_group(standard="control", exp=self.sample3)
        dfs = self.diff_freq_table(df1, df2, df3, self.sample1, self.sample2, self.sample3)
        print "Number of sites", dfs.shape[0]
        return dfs


    def diff_freq_table(self, df1=None, df2=None, df3=None, df1name="", df2name="", df3name=""):
        df1_diff = df1[["diff"]]
        df1_diff.dropna(axis=0, inplace=True)
        df2_diff = df2[["diff"]]
        df2_diff.dropna(axis=0, inplace=True)
        df3_diff = df3[["diff"]]
        df3_diff.dropna(axis=0, inplace=True)
        df12 = df1_diff.join(df2_diff, how="inner", lsuffix="_"+df1name, rsuffix="_"+df2name)
        df123 = df12.join(df3_diff, how="inner", rsuffix = "_"+df3name)
        print df123.head()
        return df123


    def significant_sites(self, cutoff1=-0.3, cutoff2=-0.3, cutoff3=-0.3):
        df = self.combine_tables()
        condition1 = df.iloc[:,0] <= cutoff1
        condition2 = df.iloc[:,1] <= cutoff2
        condition3 = df.iloc[:,2] <= cutoff3
        df_sig = df[condition1 & condition2 & condition3]
        print "number of significant sites ", df_sig.shape[0]
        return fg.add_annot(df_sig)


def main():
    sample3 = Three_sample_corr("EZH2", "ADAR1", "ADAR2")
    three_sample_sig_sites = sample3.significant_sites()
    # print three_sample_sig_sites
    pass


if __name__ == "__main__":
    main()