import Class_combine_freq_table as cmb
import Class_find_genes as fg
from matplotlib_venn import venn2, venn3
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
        print "total valued sites in all three samples", dfs.shape[0]
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
        return df123


    def significant_sites(self, cutoff1=-0.3, cutoff2=-0.3, cutoff3=-0.3): # all of them decrease <= (-0.3)
        df = self.combine_tables()
        condition1 = df.iloc[:,0] <= cutoff1
        condition2 = df.iloc[:,1] <= cutoff2
        condition3 = df.iloc[:,2] <= cutoff3
        df_sig = df[condition1 & condition2 & condition3]
        print "number of significant sites ", df_sig.shape[0]
        return fg.add_annot(df_sig)  # add annotation


    def get_site_list(self, sample, diffCutoff = -0.1):
        df = cmb.organize_group(standard="control", exp=sample)
        df_sub = df.ix[(df.ix[:, "diff"] <= diffCutoff),]
        return df_sub.index.tolist()


    def draw_venn3(self,lst1, lst2, lst3, name1, name2, name3):
        venn3([set(lst1), set(lst2), set(lst3)], set_labels=[name1, name2, name3])
        plt.show()
        pass


    def venn_3sample(self, vennCut1= (-0.1), vennCut2 = -0.1, vennCut3 = -0.1):
        """
        :param vennCut1:
        :param vennCut2:
        :param vennCut3:
        :return: Venn diagram, and each section save as a xls file
        """
        list1 = self.get_site_list(sample = self.sample1, diffCutoff = vennCut1)
        list2 = self.get_site_list(sample = self.sample2, diffCutoff = vennCut2)
        list3 = self.get_site_list(sample = self.sample3, diffCutoff = vennCut3)
        self.draw_venn3(list1, list2, list3, self.sample1, self.sample2, self.sample3)

        sample1_sample2_inter = set(list1).intersection(set(list2))
        sample2_sample3_inter = set(list2).intersection(set(list3))
        sample1_sample3_inter = set(list1).intersection(set(list3))
        all_inter = sample1_sample2_inter.intersection(set(list3))
        print "first intersection ", len(sample1_sample2_inter), self.sample1, self.sample2
        print "second intersection ", len(sample2_sample3_inter), self.sample2, self.sample3
        print "third intersection ", len(sample1_sample3_inter), self.sample1, self.sample3
        print "all intersection ", len(all_inter), self.sample1, self.sample2, self.sample3
        sample1_sample2_inter_df = self.index_df(sample1_sample2_inter)
        sample2_sample3_inter_df = self.index_df(sample2_sample3_inter)
        sample1_sample3_inter_df = self.index_df(sample1_sample3_inter)
        sample_all_inter_df = self.index_df(all_inter)
        # sample1_sample2_inter_df.to_csv("../tempResults/"+self.sample1 + "_" + self.sample2 + "_intersection.xls", sep="\t")
        # sample2_sample3_inter_df.to_csv("../tempResults/" + self.sample2 + "_" + self.sample3 + "_intersection.xls",sep="\t")
        # sample1_sample3_inter_df.to_csv("../tempResults/" + self.sample1 + "_" + self.sample3 + "_intersection.xls",sep="\t")
        # sample_all_inter_df.to_csv("../IntersectionResults/all_intersection.xls", sep="\t")
        pass


    def index_df(self, listSite):
        df = self.combine_tables()
        sub_df = df.loc[df.index.intersection(listSite)]
        sub_df = fg.add_annot(sub_df)
        return sub_df


def main():
    sample3 = Three_sample_corr("EZH2", "ADAR1", "ADAR2")
    # three_sample_sig_sites = sample3.significant_sites()
    sample3.venn_3sample()
    # print three_sample_sig_sites
    pass


if __name__ == "__main__":
    main()