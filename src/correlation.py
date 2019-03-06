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

    def combine_two_df(self):
        df1 = cmb.organize_group(standard="control", exp=self.sample1)
        df2 = cmb.organize_group(standard="control", exp=self.sample2)
        df1_df2 = self.diff_freq_table(df1, df2, self.sample1, self.sample2)
        print "Number of sites", df1_df2.shape[0]
        return df1_df2


    def diff_freq_table(self, df1=None, df2=None, df1name="", df2name="", definedSites = False):
        """
        first df remove nan, df2 replace nan to 0; join them left.
        sub function used by combine_two_df
        defineSites: True, if ADAR cutoff is used as 0.10
        """
        df1_diff = df1[["diff"]]
        df1_diff.dropna(axis=0, inplace=True)
        if definedSites:
            df2 = df2[df2["diff"] <= (-0.1)]
        df2_diff = df2[["diff"]]
        df2_diff.dropna(axis=0, inplace=True)
        df = df1_diff.join(df2_diff, how="inner", lsuffix="_"+df1name, rsuffix="_"+df2name)
        return df


    def correlation_value(self):
        df = self.combine_two_df()
        p_pearson =  stat.pearsonr(df.iloc[:,0],df.iloc[:,1])
        p_spearman = stat.spearmanr(df.iloc[:,0], df.iloc[:,1])
        print "Pearson correlation: ", p_pearson
        print "Pearson p value: %.8f" % p_pearson[1]
        print  p_spearman
        print "Spearman p value: %.8f" %p_spearman[1]
        pass


    def plot_correlation(self):
        df = self.combine_two_df()
        x = df.iloc[:,0]
        y = df.iloc[:,1]
        plt.scatter(x, y)
        plt.xlabel(df.columns[0])
        plt.ylabel(df.columns[1])
        plt.show()
        pass


    def significant_sites(self, cutoff1=-0.3, cutoff2=-0.3):
        df = self.combine_two_df()
        condition1 = df.iloc[:,0] <= cutoff1
        condition2 = df.iloc[:,1] <= cutoff2
        df_sig = df[condition1 & condition2]
        print "number of significant sites ", df_sig.shape[0]
        return fg.add_annot(df_sig)


    def two_sample_correlation(self):
        self.plot_correlation()
        self.correlation_value()
        self.significant_sites()
        pass


def main():
    # EZH2_ADAR1 = Two_sample_corr("EZH2", "ADAR1")
    # EZH2_ADAR1.two_sample_correlation()
    # EZH2_ADAR1_sigSites = EZH2_ADAR1.significant_sites(-0.4, -0.4)
    # print EZH2_ADAR1_sigSites

    # EZH2_ADAR2 = Two_sample_corr("EZH2", "ADAR2")
    # EZH2_ADAR2.two_sample_correlation()
    # EZH2_ADAR2_sigSites = EZH2_ADAR2.significant_sites(-0.4, -0.4)
    # print EZH2_ADAR2_sigSites

    ADAR1_ADAR2 = Two_sample_corr("ADAR1", "ADAR2")
    ADAR1_ADAR2.two_sample_correlation()
    # ADAR1_ADAR2_sigSites = ADAR1_ADAR2.significant_sites(-0.4, -0.4)
    # print ADAR1_ADAR2_sigSites
    pass


if __name__ == "__main__":
    main()