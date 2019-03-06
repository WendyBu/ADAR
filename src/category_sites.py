import Class_combine_freq_table as cmb
import Class_find_genes as fg
from matplotlib_venn import venn2, venn3
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np


class Category_site:
    def __init__(self, df):
        self.df = df
        pass

    def count_types(self):
        groups = self.df.groupby("annot1").size()
        print groups
        pass


    def get_3UTR(self):
        df_3UTR = (self.df)["annot1"] == "3UTR"
        df_3UTR = (self.df).ix[df_3UTR,]
        return df_3UTR


    def get_5UTR(self):
        df_5UTR = (self.df)["annot1"] == "5UTR"
        df_5UTR = (self.df).ix[df_5UTR,]
        return df_5UTR


    def get_Intron(self):
        df_intron = (self.df)["annot1"] == "intronic"
        df_intron = (self.df).ix[df_intron,]
        return df_intron


def main():
    # inputFile = "../tempResults/EZH2_ADAR1_intersection.xls"
    # df = pd.read_csv(inputFile, sep="\t", index_col=0)
    # df = df.sort_values("diff_EZH2")
    # EZH2_ADAR1 = Category_site(df)
    # EZH2_ADAR1.count_types()
    # EZH2_ADAR1_3UTR = EZH2_ADAR1.get_3UTR()
    # EZH2_ADAR1_3UTR.to_csv("../Compare_category/EZH2_ADAR1_3UTR.xls", sep="\t")
    # EZH2_ADAR1_5UTR = EZH2_ADAR1.get_5UTR()
    # EZH2_ADAR1_5UTR.to_csv("../Compare_category/EZH2_ADAR1_5UTR.xls", sep="\t")
    # EZH2_ADAR1_intron = EZH2_ADAR1.get_Intron()
    # EZH2_ADAR1_intron.to_csv("../Compare_category/EZH2_ADAR1_intron.xls", sep="\t")

    # inputFile = "../tempResults/EZH2_ADAR2_intersection.xls"
    # df = pd.read_csv(inputFile, sep="\t", index_col=0)
    # df = df.sort_values("diff_EZH2")
    # EZH2_ADAR2 = Category_site(df)
    # EZH2_ADAR2.count_types()
    # EZH2_ADAR2_3UTR = EZH2_ADAR2.get_3UTR()
    # EZH2_ADAR2_3UTR.to_csv("../Compare_category/EZH2_ADAR2_3UTR.xls", sep="\t")
    # EZH2_ADAR2_5UTR = EZH2_ADAR2.get_5UTR()
    # EZH2_ADAR2_5UTR.to_csv("../Compare_category/EZH2_ADAR2_5UTR.xls", sep="\t")
    # EZH2_ADAR2_intron = EZH2_ADAR2.get_Intron()
    # EZH2_ADAR2_intron.to_csv("../Compare_category/EZH2_ADAR2_intron.xls", sep="\t")

    # inputFile = "../tempResults/all_intersection.xls"
    # df = pd.read_csv(inputFile, sep="\t", index_col=0)
    # df = df.sort_values("diff_EZH2")
    # all_inter = Category_site(df)
    # all_inter.count_types()
    # all_inter_3UTR = all_inter.get_3UTR()
    # all_inter_3UTR.to_csv("../Compare_category/all_inter_3UTR.xls", sep="\t")
    # all_inter_5UTR = all_inter.get_5UTR()
    # all_inter_5UTR.to_csv("../Compare_category/all_inter_5UTR.xls", sep="\t")
    # all_inter_intron = all_inter.get_Intron()
    # all_inter_intron.to_csv("../Compare_category/all_inter_intron.xls", sep="\t")

    inputFile = "../tempResults/ADAR1_ADAR2_intersection.xls"
    df = pd.read_csv(inputFile, sep="\t", index_col=0)
    df = df.sort_values("diff_ADAR1")
    ADAR1_ADAR2 = Category_site(df)
    ADAR1_ADAR2.count_types()
    # ADAR1_ADAR2_3UTR = ADAR1_ADAR2.get_3UTR()
    # ADAR1_ADAR2_3UTR.to_csv("../Compare_category/ADAR1_ADAR2_3UTR.xls", sep="\t")
    # ADAR1_ADAR2_5UTR = ADAR1_ADAR2.get_5UTR()
    # ADAR1_ADAR2_5UTR.to_csv("../Compare_category/ADAR1_ADAR2_5UTR.xls", sep="\t")
    # ADAR1_ADAR2_intron = ADAR1_ADAR2.get_Intron()
    # ADAR1_ADAR2_intron.to_csv("../Compare_category/ADAR1_ADAR2_intron.xls", sep="\t")
    pass


if __name__ == "__main__":
    main()