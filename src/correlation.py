import Class_combine_freq_table as cmb
import Class_define_cutoff as cf
import Class_find_genes as fg
import numpy as np
import scipy.stats as stat
from matplotlib import pyplot as plt


class Two_sample_corr:
    def __init__(self, sample1, sample2):
        self.sample1 = sample1
        self.sample2 = sample2


    def correlation_value(self, df):
        p_pearson =  stat.pearsonr(df.iloc[:,0],df.iloc[:,1])
        p_spearman = stat.spearmanr(df.iloc[:,0], df.iloc[:,1])
        print "Pearson correlation: ", p_pearson
        print "Pearson p value: %.8f" % p_pearson[1]
        print  p_spearman
        print "Spearman p value: %.8f" %p_spearman[1]
        pass


    def plot_correlation(self, df):
        x = df.iloc[:,0]
        y = df.iloc[:,1]
        plt.scatter(x, y)
        plt.xlabel(df.columns[0])
        plt.ylabel(df.columns[1])
        plt.show()
        pass


    def diff_freq_table(self, df1=None, df2=None, df1name="", df2name=""):
        """
        first df remove nan, df2 replace nan to 0; join them left.
        """
        df1_diff = df1[["diff"]]
        df1_diff.dropna(axis=0, inplace=True)
        df2_diff = df2[["diff"]]
        df2_diff.dropna(axis=0, inplace=True)
        df = df1_diff.join(df2_diff, how="inner", lsuffix="_"+df1name, rsuffix="_"+df2name)
        return df


    def two_sample_correlation(self):
        df1 = cmb.organize_group(standard="control", exp=self.sample1)
        df2 = cmb.organize_group(standard="control", exp=self.sample2)
        df1_df2 = self.diff_freq_table(df1, df2, self.sample1, self.sample2)
        print "Number of sites", df1_df2.shape[0]
        self.plot_correlation(df1_df2)
        self.correlation_value(df1_df2)
        pass


def main():
    EZH2_ADAR1_compare = Two_sample_corr("EZH2", "ADAR1")
    EZH2_ADAR1_compare.two_sample_correlation()

    EZH2_ADAR1_compare = Two_sample_corr("EZH2", "ADAR2")
    EZH2_ADAR1_compare.two_sample_correlation()

    EZH2_ADAR1_compare = Two_sample_corr("ADAR1", "ADAR2")
    EZH2_ADAR1_compare.two_sample_correlation()
    pass


if __name__ == "__main__":
    main()