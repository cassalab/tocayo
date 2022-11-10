import pandas as pd


if __name__ == "__main__":

	pd.options.mode.chained_assignment = None

	classificationList = ["P/LP", "B/LB"]
	classificationList_otr = ["B/LB", "P/LP"]
	classificationList_v2 = ["plp", "blb"]
	classificationList_v2_otr = ["blb", "plp"]

	for i in range(len(classificationList)):

		clas = classificationList[i]
		clas_otr = classificationList_otr[i]
		clas_v2 = classificationList_v2[i]
		clas_v2_otr = classificationList_v2_otr[i]

		nameDf = pd.read_csv("/net/data/aasubs/clinvar_only/dup_files/ps1_all_subs.csv")

		df = pd.read_csv(f"/net/data/aasubs/clinvar_only/annotated/strong_annotated_dup_subs.csv")
		df = df[df["simple_name"].isin(nameDf["name"])].sort_values(by = ["vus_aa_sub_name"]).reset_index(drop = True)
		df.drop(df.filter(regex="Unname"),axis=1, inplace=True)

		ndDf = df.drop_duplicates(subset = ["aa_sub_name"]).reset_index(drop = True)

		finalClasSepDf = pd.DataFrame(columns = df.columns)
		finalVusSepDf = pd.DataFrame(columns = df.columns)

		for _, uniq_entry in ndDf.iterrows():
			
			all_entries = df[df["aa_sub_name"] == uniq_entry["aa_sub_name"]].reset_index(drop = True)

			if len(all_entries[all_entries["simple_annot"] == "VUS"]) > 0 and len(all_entries[all_entries["simple_annot"] == clas]) > 0 and len(all_entries[all_entries["simple_annot"] == clas_otr]) == 0:

				plp_plp_df = all_entries[all_entries["simple_annot"] == clas].reset_index(drop = True)
				plp_vus_df = all_entries[all_entries["simple_annot"] == "VUS"].reset_index(drop = True)

				for x, plp_entry in plp_plp_df.iterrows():
					for y, vus_entry in plp_vus_df.iterrows():
						finalClasSepDf = pd.concat([finalClasSepDf, plp_plp_df.iloc[[x]]], ignore_index = True)
						finalVusSepDf = pd.concat([finalVusSepDf, plp_vus_df.iloc[[y]]], ignore_index = True)

		finalClasSepDf.drop(finalClasSepDf.filter(regex="Unname"),axis=1, inplace=True)
		finalVusSepDf.drop(finalVusSepDf.filter(regex="Unname"),axis=1, inplace=True)

		clas_columns_sep = []
		for column in list(finalClasSepDf.columns):
			clas_columns_sep.append(clas_v2 + "_" + str(column))
		finalClasSepDf.columns = clas_columns_sep

		vus_columns_sep = []
		for column in list(finalVusSepDf.columns):
			vus_columns_sep.append("vus_" + str(column))
		finalVusSepDf.columns = vus_columns_sep

		finalSepDifDf = pd.concat([finalClasSepDf, finalVusSepDf], axis = 1)
		finalSepDifDf = finalSepDifDf.sort_values(by = ["vus_aa_sub_name"]).reset_index(drop = True)

		finalSepDifDf.to_csv(f"/net/data/aasubs/clinvar_only/results/strong_{clas_v2}_analysis.csv")

