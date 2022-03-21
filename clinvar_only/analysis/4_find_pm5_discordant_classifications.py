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

		nameDf = pd.read_csv("/net/data/aasubs/clinvar_only/dup_files/pm5_all_subs.csv")

		df = pd.read_csv(f"/net/data/aasubs/clinvar_only/annotated/moderate_annotated_dup_subs.csv")
		df = df[df["simple_name"].isin(nameDf["name"])].sort_values(by = ["simple_name"]).reset_index(drop = True)
		df.drop(df.filter(regex="Unname"),axis=1, inplace=True)

		ndDf = df.drop_duplicates(subset = ["aa_half_name"])

		finalClasSepDf = pd.DataFrame(columns = df.columns)
		finalVusSepDf = pd.DataFrame(columns = df.columns)

		for y, uniq_entry in ndDf.iterrows():

			if len(all_entries[all_entries["simple_annot"] == "VUS"]) > 0 and len(all_entries[all_entries["simple_annot"] == clas]) > 0 and len(all_entries[all_entries["simple_annot"] == clas_otr]) == 0:

				plp_plp_df = all_entries[all_entries["simple_annot"] == clas].reset_index(drop = True)
				plp_vus_df = all_entries[all_entries["simple_annot"] == "VUS"].reset_index(drop = True)

				for x, plp_entry in plp_plp_df.iterrows():
					for y, vus_entry in plp_vus_df.iterrows():
						if plp_entry["aa_sub_name"] != vus_entry["aa_sub_name"]:
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
		finalSepDifDfND = finalSepDifDf.drop_duplicates(subset = [f"{clas_v2}_aa_sub_name"]).drop_duplicates(subset = ["vus_aa_sub_name"])

		finalSepDifDf.to_csv(f"/net/data/aasubs/clinvar_only/results_new/moderate_{clas_v2}_analysis.csv")

