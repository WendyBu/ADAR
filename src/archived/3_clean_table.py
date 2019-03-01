import pandas as pd
import glob
import os.path
pd.set_option("display.max_column", 100)

# keep all the items if stand is + and transition is AG
# keep all the items if strand is - and transition is TC
for file in glob.glob("../AGTC_ref/*"):
    df = pd.read_csv(file, sep="\t")
    negative_strand = (df.strand=="-") & (df.AllSubs == "TC")
    positive_strand = (df.strand=="+") & (df.AllSubs == "AG")
    df_correct_strand = df[negative_strand|positive_strand]

    df_clean = df_correct_strand.loc[:, ["ID", "Region", "Position", "BaseCount[A,C,G,T]", "AllSubs", "Frequency",\
                                          "gene", "strand", "identifier", "annot1", "annot2" , "alu?", "non_alu_repetitive?"]]
    df_clean.to_csv("../clean_AGTC/" + os.path.basename(file), sep="\t")

