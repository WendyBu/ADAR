"""
add identifier to the download editing reference
reference downloaded from http://rnaedit.com/download/
hg19 human
"""

import pandas as pd
pd.set_option("display.max_column", 100)

data = pd.read_csv("../database/Human_AG_all_hg19_v2.txt", sep="\t")
data["identifier"] = data['chromosome'] + "_" + data["position"].astype(str)
data.set_index("identifier", inplace=True)
data.to_csv("../database/hg19_AG_editing_reference.txt", sep="\t")