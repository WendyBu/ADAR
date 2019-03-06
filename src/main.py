import Class_combine_freq_table as cmb
import Class_define_cutoff as cf
import Class_find_genes as fg


def main():
    # df = cmb.organize_group(standard="control", exp="ADAR2")  # get the frequency table
    # standard_cutoff, siteList, length = cf.ROC(df)            #  find the cutoff using the exp label
    # print "final standard cutoff ", standard_cutoff
    # print "number of adar editing sites ", length
    df_target = cmb.organize_group(standard="control", exp="EZH2")   # get the frequency table
    # sample_table = fg.find_all_genes(df_target, siteList)            # get the gene list fit all three sample labels
    # sample_table.to_csv("../results/control_ezh2_adar.xls", sep="\t")
    pass


if __name__ == "__main__":
    main()