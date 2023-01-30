import pandas as pd

evidenceList = ["strong", "moderate"]
classificationList = ["P/LP", "B/LB"]
classificationList_v2 = ["plp", "blb"]

for evidence in evidenceList:
	for i in range(len(classificationList)):

		clas = classificationList[i]
		clas_v2 = classificationList_v2[i]
		
		df = pd.read_csv(f"/Users/vineel/Desktop/dbnsfp_{evidence}_{clas_v2}_analysis.csv")
		df["cv_simple_name"] = df["cv_Chromosome"].map(str) + "-" + df["cv_Start"].map(str) + "-" + df["cv_AlternateAlleleVCF"].map(str)
		df["db_simple_name"] = df["db_chr"].map(str) + "-" + df["db_pos"].map(str) + "-" + df["db_alt"].map(str)
		
		#SpliceAI data retrieved from VEP Web Interface
		addCvDf = pd.read_csv(f"/Users/vineel/Documents/dbnsfp_sa_vep/dbnsfp_{evidence}_{clas_v2}_cv_output.txt", sep = "\t")
		addDbDf = pd.read_csv(f"/Users/vineel/Documents/dbnsfp_sa_vep/dbnsfp_{evidence}_{clas_v2}_db_output.txt", sep = "\t")

		addCvDf["vep_chrom"] = addCvDf["Location"].str.split(":").str.get(0)
		addCvDf["vep_pos"] = addCvDf["Location"].str.split("-").str.get(1)

		addDbDf["vep_chrom"] = addDbDf["Location"].str.split(":").str.get(0)
		addDbDf["vep_pos"] = addDbDf["Location"].str.split("-").str.get(1)

		addCvDf["simple_name"] = addCvDf["vep_chrom"].map(str) + "-" + addCvDf["vep_pos"].map(str) + "-" + addCvDf["Allele"].map(str)
		addDbDf["simple_name"] = addDbDf["vep_chrom"].map(str) + "-" + addDbDf["vep_pos"].map(str) + "-" + addDbDf["Allele"].map(str)

		addCvDf = addCvDf[(addCvDf["SpliceAI_pred_DS_AG"] != "-") | (addCvDf["SpliceAI_pred_DS_AL"] != "-") | (addCvDf["SpliceAI_pred_DS_DG"] != "-") | (addCvDf["SpliceAI_pred_DS_DL"] != "-")]
		addCvDf = addCvDf.drop_duplicates(subset = ["simple_name"])

		addDbDf = addDbDf[(addDbDf["SpliceAI_pred_DS_AG"] != "-") | (addDbDf["SpliceAI_pred_DS_AL"] != "-") | (addDbDf["SpliceAI_pred_DS_DG"] != "-") | (addDbDf["SpliceAI_pred_DS_DL"] != "-")]
		addDbDf = addDbDf.drop_duplicates(subset = ["simple_name"])

		addCvDf["sa_score"] = addCvDf[["SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL"]].astype(float).max(axis = 1)
		addDbDf["sa_score"] = addDbDf[["SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL"]].astype(float).max(axis = 1)

		df = df[(df["cv_simple_name"].isin(addCvDf["simple_name"])) & (df["db_simple_name"].isin(addDbDf["simple_name"]))]

		addCvDf = addCvDf[addCvDf["simple_name"].isin(df["cv_simple_name"])].sort_values(by = ["simple_name"]).reset_index(drop = True)
		addDbDf = addDbDf[addDbDf["simple_name"].isin(df["db_simple_name"])].sort_values(by = ["simple_name"]).reset_index(drop = True)

		df = df.sort_values(by = ["cv_simple_name"]).reset_index(drop = True)
		cvCI = df.drop_duplicates(subset = ["cv_simple_name"])
		addCvDf = addCvDf.set_index(cvCI.index)
		addCvDf = addCvDf.reindex(df.index, method = "ffill").sort_values(by = ["simple_name"]).reset_index(drop = True)
		df["cv_sa_score"] = addCvDf["sa_score"]

		df = df.sort_values(by = ["db_simple_name"]).reset_index(drop = True)
		dbCI = df.drop_duplicates(subset = ["db_simple_name"])
		addDbDf = addDbDf.set_index(dbCI.index)
		addDbDf = addDbDf.reindex(df.index, method = "ffill").sort_values(by = ["simple_name"]).reset_index(drop = True)
		df["db_sa_score"] = addDbDf["sa_score"]

		df["cv_simple_name"] = df["cv_Chromosome"].map(str) + "-" + df["cv_Start"].map(str) + "-" + df["cv_ReferenceAlleleVCF"].map(str) + "-" + df["cv_AlternateAlleleVCF"].map(str)
		df["db_simple_name"] = df["db_chr"].map(str) + "-" + df["db_pos"].map(str) + "-" + df["db_ref"].map(str) + "-" + df["db_alt"].map(str)
		
		df.drop(df.filter(regex="Unname"),axis=1, inplace=True)
		df = df.sort_values(by = ["cv_aa_sub_name"]).reset_index(drop = True)
		df.to_csv(f"/Users/vineel/Documents/splicing/dbnsfp/dbnsfp_{evidence}_{clas_v2}_spliceai.csv", index = False)

