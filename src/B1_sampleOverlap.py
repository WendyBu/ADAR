import pandas as pd
import glob
import os.path
from matplotlib import pyplot as plt
import numpy as np
from matplotlib_venn import venn2
from matplotlib_venn import venn3, venn3_circles
pd.set_option("display.max_column", 100)





def diff_genes(file1, file2):
    df1 = pd.read_csv(file1, sep="\t")
    df2 = pd.read_csv(file2, sep="\t")
    df1_list = df1.gene.tolist()
    df2_list = df2.gene.tolist()
    # venn3([set(control), set(EZH2sh1), set(EZH2sh2)], set_labels=["control", "EZH2sh1", "EZH2sh2"])
    # venn3([set(control), set(ADARsh1), set(ADARsh2)], set_labels=["control", "ADARsh1", "ADARsh2"])
    venn2([set(df1_list), set(df2_list)], set_labels=[os.path.basename(file1).split("_")[0], os.path.basename(file2).split("_")[0]])
    plt.show()
    pass


def compare_frequence_duplicates(file1, file2):
    df1 = pd.read_csv(file1, sep="\t")
    df2 = pd.read_csv(file2, sep="\t")
    df1.set_index("ID", inplace=True)
    df2.set_index("ID", inplace=True)
    df1 = df1.loc[:, ["Frequency", "BaseCount[A,C,G,T]"]]
    df2 = df2.loc[:, ["Frequency", "BaseCount[A,C,G,T]"]]
    compare_2df = df1.join(df2, how="outer", rsuffix="_2ndsample")
    print compare_2df.head(100)
    pass


def generate_list_sample(file1, file2):
    df1 = pd.read_csv(file1, sep="\t")
    df2 = pd.read_csv(file2, sep="\t")
    df1_list = df1.gene.tolist()
    print "df1 length", len(df1_list)
    df2_list = df2.gene.tolist()
    print "df2 length", len(df2_list)
    df = df1_list + df2_list
    print "df length", len(df)
    df_set = set(df)
    print "df set length", len(df_set)
    return df_set



def main():
    # check sample duplicates difference
    # diff_genes("../clean_AGTC/A1_AGTC.txt", "../clean_AGTC/A2_AGTC.txt")
    # diff_genes("../clean_AGTC/A3_AGTC.txt", "../clean_AGTC/A4_AGTC.txt")
    # diff_genes("../clean_AGTC/A5_AGTC.txt", "../clean_AGTC/A6_AGTC.txt")
    # diff_genes("../clean_AGTC/A7_AGTC.txt", "../clean_AGTC/A8_AGTC.txt")
    # diff_genes("../clean_AGTC/A9_AGTC.txt", "../clean_AGTC/A10_AGTC.txt")
    # diff_genes("../clean_AGTC/A11_AGTC.txt", "../clean_AGTC/A12_AGTC.txt")
    # compare_frequence_duplicates("../clean_AGTC/A1_AGTC.txt", "../clean_AGTC/A2_AGTC.txt")
    control_unique_site = generate_list_sample("../clean_AGTC/A1_AGTC.txt", "../clean_AGTC/A2_AGTC.txt")
    EZH2sh1_unique_site = generate_list_sample("../clean_AGTC/A3_AGTC.txt", "../clean_AGTC/A4_AGTC.txt")
    EZH2sh2_unique_site = generate_list_sample("../clean_AGTC/A5_AGTC.txt", "../clean_AGTC/A6_AGTC.txt")
    ADARsh1_unique_site = generate_list_sample("../clean_AGTC/A7_AGTC.txt", "../clean_AGTC/A8_AGTC.txt")
    ADARsh2_unique_site = generate_list_sample("../clean_AGTC/A9_AGTC.txt", "../clean_AGTC/A10_AGTC.txt")
    ADAR2sh_unique_site = generate_list_sample("../clean_AGTC/A11_AGTC.txt", "../clean_AGTC/A12_AGTC.txt")
    # venn2([control_unique_site, ADARsh1_unique_site], set_labels=["control", "ADARsh1"])
    #venn3([control_unique_site, ADARsh1_unique_site, ADAR2sh_unique_site], set_labels=["control", "ADARsh1", "ADARsh2"])
    venn3([control_unique_site, EZH2sh1_unique_site, EZH2sh2_unique_site], set_labels=["control", "EZH2sh1", "EZH2sh2"])
    # venn2([control_unique_site, ADAR2sh_unique_site], set_labels=["control", "ADAR2sh"])
    plt.show()





    pass


if __name__ == "__main__":
    main()