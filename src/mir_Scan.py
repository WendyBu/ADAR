import pandas as pd
import glob
pd.set_option('display.max_columns', 100)
import os.path


def get_ref(file):
    df = pd.read_csv(file, sep="\t", header=None, low_memory=False)
    df2 = pd.DataFrame()
    df2[["chr", "start", "end", "paired"]] = df.iloc[:, 0:4]
    return df2


def findMatchbyLine(position, ref_df):
    return ref_df[(ref_df["start"]<position) & (ref_df["end"]>position)]


def matchINDEX(df, ref):
    Result = pd.DataFrame()
    for index, row in df.iterrows():
        matchedRefRow = findMatchbyLine(row["pos"], ref)
        if not matchedRefRow.empty:
            matchedRefRow["EditingPosition"] = index
            Result = Result.append(matchedRefRow)
    return Result


def match(refFile, sample):
    ref = get_ref(refFile)
    sample_miR = matchINDEX(sample, ref)
    return sample_miR


def add_frequency(df_3UTR, df):
    df_result = df_3UTR.merge(df[["diff_EZH2", "diff_ADAR1", "diff"]], left_on="EditingPosition", right_index=True, how="left")
    return df_result


def main():
    for sampleFile in glob.glob("../Compare_category/*.xls"):
        sample = pd.read_csv(sampleFile, sep="\t", index_col=0)
        sample[["chr", "pos"]] = sample.index.to_series().str.split("_", expand=True)
        sample["pos"] = pd.to_numeric(sample["pos"])
        results = pd.DataFrame()
        for refFile in glob.glob("../miR_all_location/*.bed"):
            miR_3UTR = match(refFile, sample)
            results = results.append(miR_3UTR)

        results = add_frequency(results, sample)
        baseFileName = os.path.basename(sampleFile)
        outputFileName = os.path.join("../sample_3UTR_sites", baseFileName )
        results.to_csv(outputFileName, sep="\t", index=False)
    pass


if __name__ == "__main__":
    main()

