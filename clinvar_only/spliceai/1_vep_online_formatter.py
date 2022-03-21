import pandas as pd
from natsort import natsorted

evidenceList = ["strong", "moderate"]
classificationList_v2 = ["plp", "blb"]

for evidence in evidenceList:
	for i in range(len(classificationList_v2)):

		clas_v2 = classificationList_v2[i]

		df = pd.read_csv(f"/Users/vineel/Desktop/{evidence}_{clas_v2}_analysis.csv")
		df[f"{clas_v2}_simple_name"] = df[f"{clas_v2}_Chromosome"].map(str) + " " + df[f"{clas_v2}_Start"].map(str) + " . " + df[f"{clas_v2}_ReferenceAlleleVCF"].map(str) + " " + df[f"{clas_v2}_AlternateAlleleVCF"].map(str) + " . . ."
		df["vus_simple_name"] = df["vus_Chromosome"].map(str) + " " + df["vus_Start"].map(str) + " . " + df["vus_ReferenceAlleleVCF"].map(str) + " " + df["vus_AlternateAlleleVCF"].map(str) + " . . ."

		clasDf = df[[f"{clas_v2}_simple_name"]]
		clasDf.columns = [0]

		vusDf = df[["vus_simple_name"]]
		vusDf.columns = [0]

		clasDf[0] = natsorted(list(clasDf[0]))
		vusDf[0] = natsorted(list(vusDf[0]))

		clasDf.to_csv(f"/Users/vineel/Documents/clinvar_sa_vep/{evidence}_{clas_v2}_{clas_v2}_input.vcf", index = False, header = None)
		vusDf.to_csv(f"/Users/vineel/Documents/clinvar_sa_vep/{evidence}_{clas_v2}_vus_input.vcf", index = False, header = None)
