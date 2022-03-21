import pandas as pd
from natsort import natsorted

evidenceList = ["strong", "moderate"]
classificationList_v2 = ["plp", "blb"]

for evidence in evidenceList:
	for i in range(len(classificationList_v2)):

		clas_v2 = classificationList_v2[i]

		df = pd.read_csv(f"/Users/vineel/Desktop/dbnsfp_{evidence}_{clas_v2}_analysis.csv")
		df["cv_simple_name"] = df["cv_Chromosome"].map(str) + " " + df["cv_Start"].map(str) + " . " + df["cv_ReferenceAlleleVCF"].map(str) + " " + df["cv_AlternateAlleleVCF"].map(str) + " . . ."
		df["db_simple_name"] = df["db_chr"].map(str) + " " + df["db_pos"].map(str) + " . " + df["db_ref"].map(str) + " " + df["db_alt"].map(str) + " . . ."

		cvDf = df[["cv_simple_name"]]
		cvDf.columns = [0]

		dbDf = df[["db_simple_name"]]
		dbDf.columns = [0]

		cvDf[0] = natsorted(list(cvDf[0]))
		dbDf[0] = natsorted(list(dbDf[0]))

		cvDf.to_csv(f"/Users/vineel/Documents/dbnsfp_sa_vep/dbnsfp_{evidence}_{clas_v2}_cv_input.vcf", index = False, header = None)
		dbDf.to_csv(f"/Users/vineel/Documents/dbnsfp_sa_vep/dbnsfp_{evidence}_{clas_v2}_db_input.vcf", index = False, header = None)
