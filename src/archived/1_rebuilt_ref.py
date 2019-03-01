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


"""
input:
chromosome	position	gene	strand	annot1	annot2	alu?	non_alu_repetitive?	conservation_chimp	conservation_rhesus	conservation_mouse
chr1	206256301	C1orf186	-	intronic	intronic	no	no	N	N	N
chr6	116991832	intergenic	-	intergenic	intergenic	no	no	N	N	N
chr7	30504355	NOD1	-	intronic	intronic	no	no	N	N	N
chr1	85127959	SSX2IP	-	Syn	Gln->Gln	no	no	N	N	N
chr15	100203261	MEF2A	+	intronic	intronic	no	no	N	N	N


Output:
identifier	chromosome	position	gene	strand	annot1	annot2	alu?	non_alu_repetitive?	conservation_chimp	conservation_rhesus	conservation_mouse
chr1_206256301	chr1	206256301	C1orf186	-	intronic	intronic	no	no	N	N	N
chr6_116991832	chr6	116991832	intergenic	-	intergenic	intergenic	no	no	N	N	N
chr7_30504355	chr7	30504355	NOD1	-	intronic	intronic	no	no	N	N	N
chr1_85127959	chr1	85127959	SSX2IP	-	Syn	Gln->Gln	no	no	N	N	N
"""

