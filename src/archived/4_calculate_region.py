import pandas as pd
import glob
import os.path
from matplotlib import pyplot as plt
import numpy as np
from matplotlib_venn import venn2
from matplotlib_venn import venn3, venn3_circles
pd.set_option("display.max_column", 100)


def count_Positive_negative_strand():
    # AG_editing_num = {}
    df_editing_num = pd.DataFrame()
    for file in glob.glob("../result3_AG_clean_table/*"):
        sampleName = os.path.basename(file).split(".")[0]
        df = pd.read_csv(file, sep="\t")
        numberEditing = df.shape[0]
        # AG_editing_num.update({sampleName:numberEditing})
        df_editing_num.loc[sampleName,"AG_number"] = numberEditing

    for file in glob.glob("../result3_TC_clean_table/*"):
        sampleName = os.path.basename(file).split(".")[0]
        df = pd.read_csv(file, sep="\t")
        numberEditing = df.shape[0]
        # AG_editing_num.update({sampleName:numberEditing})
        df_editing_num.loc[sampleName, "TC_number"] = numberEditing
    df_editing_num.to_csv("../result4_Stats/AG_TC_number.csv", sep="\t")
    pass


def region_counts():
    # df_region = pd.DataFrame()
    for file in glob.glob("../result3_TC_clean_table/*"):
        sampleName = os.path.basename(file).split(".")[0]
        df = pd.read_csv(file, sep="\t")
        print  "group "+ sampleName, df.groupby("annot1").count()
    pass


def get_non_sys():
    for file in glob.glob("../result3_AG_clean_table/*"):
        sampleName = os.path.basename(file).split(".")[0]
        df = pd.read_csv(file, sep="\t")
        df_nonSys = df[df.annot1 == "Nonsyn"]
        df_nonSys.sort_values(by="gene", ascending=True, inplace=True)
        df_nonSys.to_csv("../NonSys/"+sampleName+"_AG", sep="\t")

    for file in glob.glob("../result3_TC_clean_table/*"):
        sampleName = os.path.basename(file).split(".")[0]
        df = pd.read_csv(file, sep="\t")
        df_nonSys = df[df.annot1 == "Nonsyn"]
        df_nonSys.sort_values(by="gene", ascending=True, inplace=True)
        df_nonSys.to_csv("../NonSys/" + sampleName + "_TC", sep="\t")
    pass


def analyze_nonSys():
    frequency_total = pd.DataFrame()
    for file in glob.glob("../NonSys/combine/*"):
        sampleName = os.path.basename(file).split(".")[0]
        df = pd.read_csv(file, sep="\t")
        df_frequency = df.loc[:, ["ID", "Frequency"]]
        df_frequency.set_index("ID", inplace=True)
        frequency_total = frequency_total.join(df_frequency, how="outer",   rsuffix = sampleName)
        frequency_total.fillna(0, inplace=True)
    # frequency_total.to_csv("../stat_nonSys/nonSys_stat.xls", sep="\t")
    pass



def analyze_nonSys2():
    for file in glob.glob("../NonSys/combine/*"):
        sampleName = os.path.basename(file).split(".")[0]
        df = pd.read_csv(file, sep="\t")
        df_frequency = df.ID.tolist()
        print sampleName, len(df_frequency)
    pass


def analyze_nonSys_gene():
    editGene = pd.DataFrame()
    for file in glob.glob("../NonSys/combine/*"):
        sampleName = os.path.basename(file).split(".")[0]
        df = pd.read_csv(file, sep="\t")
        geneList = df.gene.tolist()
        with open("../NonSys/geneEdited.txt", "a") as f:
            f.write("%s\t" % sampleName)
            f.write("%s\n" % geneList)

    pass

def analyze_nonSys_merge_gene(input1, input2, sampleName):
    df1 = pd.read_csv(input1, sep= "\t")
    sample1 = df1.gene
    sample1.drop_duplicates(inplace=True)
    df2 = pd.read_csv(input2, sep="\t")
    sample2 = df2.gene
    sample2.drop_duplicates(inplace=True)
    df = pd.concat([sample1, sample2], axis=1)
    df.to_csv("../NonSys/new/" + sampleName + ".txt", sep="\t", index=None)
    df1_l = df1.gene.tolist()
    df2_l = df2.gene.tolist()
    mergeL = list(set(df1_l).union(set(df2_l)))
    print len(mergeL)
    print mergeL
    outputFile = "../NonSys/new/"+sampleName + "_combine.txt"
    with open(outputFile, "a")as f:
        for item in mergeL:
            f.write("%s\n" % item)
    pass


def diff_genes():
    df = pd.read_csv("../NonSys/new/All_editing_genes_summary.txt", sep="\t")
    print df.head()
    control = df.Control.tolist()
    EZH2sh1 = df.EZH2sh1.tolist()
    EZH2sh2 = df.EZH2sh2.tolist()
    ADARsh1 = df.ADAR1sh1.tolist()
    ADARsh2 = df.ADAR1sh2.tolist()
    ADAR2 = df.ADAR2.tolist()
    # venn3([set(control), set(EZH2sh1), set(EZH2sh2)], set_labels=["control", "EZH2sh1", "EZH2sh2"])
    # venn3([set(control), set(ADARsh1), set(ADARsh2)], set_labels=["control", "ADARsh1", "ADARsh2"])
    venn2([set(control), set(ADAR2)], set_labels=["control", "ADAR2"])
    plt.show()
    pass



def main():
    # step1
    # count_Positive_negative_strand()
    # step2
    # region_counts()
    # get the non-sys
    # get_non_sys()
    # analyze_nonSys()

    # input1 = "../NonSys/combine/A11.txt"
    # input2 = "../NonSys/combine/A12.txt"
    # sampleName = "ADAR2"
    # analyze_nonSys_merge_gene(input1, input2, sampleName)

    diff_genes()
    pass


if __name__ == "__main__":
    main()