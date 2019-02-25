import pandas as pd
import glob
import os.path
from matplotlib import pyplot as plt
import numpy as np
from matplotlib_venn import venn2
from matplotlib_venn import venn3, venn3_circles
pd.set_option("display.max_column", 100)


def average_frequency():
    for file in glob.glob("../clean_AGTC/*"):
        df = pd.read_csv(file, sep="\t")
        # print df[["Frequency"]]
        freq = df.loc[:, "Frequency"].mean()
        print file
        print freq
    pass


def exon_average_frequency():
    for file in glob.glob("../clean_AGTC/*"):
        df = pd.read_csv(file, sep="\t")
        # print df.annot1.unique()
        df_intron = df[df.annot1 != "intergenic"]
        freq_intron =  df_intron.loc[:, "Frequency"].mean()

        print file
        print freq_intron
    pass




def main():
    # average_frequency()
    # exon_average_frequency()



    pass


if __name__ == "__main__":
    main()