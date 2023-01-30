import pandas as pd

evidenceList = ["strong", "moderate"]
classificationList = ["P/LP", "B/LB"]
classificationList_v2 = ["plp", "blb"]

for evidence in evidenceList:
	for i in range(len(classificationList)):

		clas = classificationList[i]
		clas_v2 = classificationList_v2[i]
    
		df = pd.read_csv(f"/Users/vineel/Desktop/{evidence}_{clas_v2}_analysis.csv")
		df[f"{clas_v2}_simple_name"] = df[f"{clas_v2}_Chromosome"].map(str) + "-" + df[f"{clas_v2}_Start"].map(str) + "-" + df[f"{clas_v2}_AlternateAlleleVCF"].map(str)
		df["vus_simple_name"] = df["vus_Chromosome"].map(str) + "-" + df["vus_Start"].map(str) + "-" + df["vus_AlternateAlleleVCF"].map(str)
		
		#SpliceAI data retrieved from VEP Web Interface
		addClasDf = pd.read_csv(f"/Users/vineel/Documents/sa_vep/{evidence}_{clas_v2}_{clas_v2}_output.txt", sep = "\t")
		addVusDf = pd.read_csv(f"/Users/vineel/Documents/sa_vep/{evidence}_{clas_v2}_vus_output.txt", sep = "\t")

		addClasDf["vep_chrom"] = addClasDf["Location"].str.split(":").str.get(0)
		addClasDf["vep_pos"] = addClasDf["Location"].str.split("-").str.get(1)

		addVusDf["vep_chrom"] = addVusDf["Location"].str.split(":").str.get(0)
		addVusDf["vep_pos"] = addVusDf["Location"].str.split("-").str.get(1)

		addClasDf["simple_name"] = addClasDf["vep_chrom"].map(str) + "-" + addClasDf["vep_pos"].map(str) + "-" + addClasDf["Allele"].map(str)
		addVusDf["simple_name"] = addVusDf["vep_chrom"].map(str) + "-" + addVusDf["vep_pos"].map(str) + "-" + addVusDf["Allele"].map(str)

		addClasDf = addClasDf[(addClasDf["SpliceAI_pred_DS_AG"] != "-") | (addClasDf["SpliceAI_pred_DS_AL"] != "-") | (addClasDf["SpliceAI_pred_DS_DG"] != "-") | (addClasDf["SpliceAI_pred_DS_DL"] != "-")]
    		addClasDf = addClasDf.drop_duplicates(subset = ["simple_name"])

		addVusDf = addVusDf[(addVusDf["SpliceAI_pred_DS_AG"] != "-") | (addVusDf["SpliceAI_pred_DS_AL"] != "-") | (addVusDf["SpliceAI_pred_DS_DG"] != "-") | (addVusDf["SpliceAI_pred_DS_DL"] != "-")]
    		addVusDf = addVusDf.drop_duplicates(subset = ["simple_name"])

		addClasDf["sa_score"] = addClasDf[["SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL"]].astype(float).max(axis = 1)
		addVusDf["sa_score"] = addVusDf[["SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL"]].astype(float).max(axis = 1)

		df = df[(df[f"{clas_v2}_simple_name"].isin(addClasDf["simple_name"])) & (df["vus_simple_name"].isin(addVusDf["simple_name"]))]

		addClasDf = addClasDf[addClasDf["simple_name"].isin(df[f"{clas_v2}_simple_name"])].sort_values(by = ["simple_name"]).reset_index(drop = True)
		addVusDf = addVusDf[addVusDf["simple_name"].isin(df["vus_simple_name"])].sort_values(by = ["simple_name"]).reset_index(drop = True)

		df = df.sort_values(by = [f"{clas_v2}_simple_name"]).reset_index(drop = True)
		clasCI = df.drop_duplicates(subset = [f"{clas_v2}_simple_name"])
		addClasDf = addClasDf.set_index(clasCI.index)
		addClasDf = addClasDf.reindex(df.index, method = "ffill").sort_values(by = ["simple_name"]).reset_index(drop = True)
		df[f"{clas_v2}_sa_score"] = addClasDf["sa_score"]

		df = df.sort_values(by = ["vus_simple_name"]).reset_index(drop = True)
		vusCI = df.drop_duplicates(subset = ["vus_simple_name"])
		addVusDf = addVusDf.set_index(vusCI.index)
		addVusDf = addVusDf.reindex(df.index, method = "ffill").sort_values(by = ["simple_name"]).reset_index(drop = True)
		df["vus_sa_score"] = addVusDf["sa_score"]

		df = df.sort_values(by = ["vus_aa_sub_name"]).reset_index(drop = True)
		df.to_csv(f"/Users/vineel/Documents/splicing/clinvar/{evidence}_{clas_v2}_spliceai.csv")

