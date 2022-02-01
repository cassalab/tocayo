import pandas as pd
import sys


def indices(lst, element):
    result = []
    offset = -1
    while True:
        try:
            offset = lst.index(element, offset+1)
        except ValueError:
            return result
        result.append(offset)


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

		df = pd.read_csv(f"/net/data/aasubs/clinvar_only/annotated/strong_annotated_dup_subs.csv")
		df = df[(df["impact"] == "MODERATE") | (df["impact"] == "HIGH")]
		df.drop(df.filter(regex="Unname"),axis=1, inplace=True)

		clasDf = df[df["simple_annot"] == clas].reset_index(drop = True)
		vusDf = df[df["simple_annot"] == "VUS"].reset_index(drop = True)
		vusSubList = list(vusDf["aa_sub_name"])
		otrDf = df[df["simple_annot"] == clas_otr].reset_index(drop = True)
		otrSubList = list(otrDf["aa_sub_name"])

		ndDf = clasDf.drop_duplicates(subset = ["aa_sub_name"])
		smDf = clasDf[(clasDf["aa_sub_name"].isin(ndDf["aa_sub_name"])) & (clasDf["simple_name"].isin(ndDf["simple_name"]) == False)]

		finalClasSepDf = pd.DataFrame(columns = df.columns)
		finalVusSepDf = pd.DataFrame(columns = df.columns)

		finalClas2SepDf = pd.DataFrame(columns = df.columns)
		finalSmSepDf = pd.DataFrame(columns = df.columns)

		finalClas3SepDf = pd.DataFrame(columns = df.columns)
		finalOtrSepDf = pd.DataFrame(columns = df.columns)

		for _, uniq_clas_entry in ndDf.iterrows():

			for y in indices(list(vusDf["aa_sub_name"]), uniq_clas_entry["aa_sub_name"]):

				vus_entry = vusDf.iloc[[y]]
				all_clas_entries = clasDf[clasDf["aa_sub_name"] == uniq_clas_entry["aa_sub_name"]].reset_index(drop = True)

				for z, clas_entry in all_clas_entries.iterrows():

					if (
						str(clas_entry["full_aa"]) == vus_entry["full_aa"].to_string(index = False, header = False)[1:] and 
						int(clas_entry["protein_position"]) == int(vus_entry["protein_position"])
						):

						clas_entry["codonpos"] = [i for i, c in enumerate(clas_entry["ref_codons"].split("/")[0]) if c.isupper()][0] + 1

						if str(clas_entry["strand"]) == "+":

							if ((
								float(clas_entry["codonpos"]) == 1 and 
								int(vus_entry["Start"]) >= int(clas_entry["Start"]) and 
								int(vus_entry["Start"]) <= int(clas_entry["Start"])+2
								)
								or 
								(
								float(clas_entry["codonpos"]) == 2 and 
								int(vus_entry["Start"]) >= int(clas_entry["Start"])-1 and 
								int(vus_entry["Start"]) <= int(clas_entry["Start"])+1
								)
								or
								(
								float(clas_entry["codonpos"]) == 3 and 
								int(vus_entry["Start"]) >= int(clas_entry["Start"])-2 and 
								int(vus_entry["Start"]) <= int(clas_entry["Start"])
								)):

								finalClasSepDf = pd.concat([finalClasSepDf, all_clas_entries.iloc[[z]]], ignore_index = True)
								finalVusSepDf = pd.concat([finalVusSepDf, vus_entry], ignore_index = True)
								
								continue

						else:

							if ((
								4-float(clas_entry["codonpos"]) == 1 and 
								int(vus_entry["Start"]) >= int(clas_entry["Start"]) and 
								int(vus_entry["Start"]) <= int(clas_entry["Start"])+2
								)
								or 
								(
								4-float(clas_entry["codonpos"]) == 2 and 
								int(vus_entry["Start"]) >= int(clas_entry["Start"])-1 and 
								int(vus_entry["Start"]) <= int(clas_entry["Start"])+1
								)
								or
								(
								4-float(clas_entry["codonpos"]) == 3 and 
								int(vus_entry["Start"]) >= int(clas_entry["Start"])-2 and 
								int(vus_entry["Start"]) <= int(clas_entry["Start"])
								)):

								finalClasSepDf = pd.concat([finalClasSepDf, all_clas_entries.iloc[[z]]], ignore_index = True)
								finalVusSepDf = pd.concat([finalVusSepDf, vus_entry], ignore_index = True)

								continue

						if str(clas_entry["vep_refseq"]) == vus_entry["vep_refseq"].to_string(index = False, header = False)[1:]:

							finalClasSepDf = pd.concat([finalClasSepDf, all_clas_entries.iloc[[z]]], ignore_index = True)
							finalVusSepDf = pd.concat([finalVusSepDf, vus_entry], ignore_index = True)

							continue

						if vus_entry["aa_sub_name"].to_string(index = False, header = False)[1:] in vusSubList:
							vusSubList.remove(vus_entry["aa_sub_name"].to_string(index = False, header = False)[1:])


			for y in indices(list(otrDf["aa_sub_name"]), uniq_clas_entry["aa_sub_name"]):

				otr_entry = otrDf.iloc[[y]]
				all_clas_entries_3 = clasDf[clasDf["aa_sub_name"] == uniq_clas_entry["aa_sub_name"]].reset_index(drop = True)

				for z, clas_entry_3 in all_clas_entries_3.iterrows():

					if (
						str(clas_entry_3["full_aa"]) == otr_entry["full_aa"].to_string(index = False, header = False)[1:] and 
						int(clas_entry_3["protein_position"]) == int(otr_entry["protein_position"])
						):

						clas_entry_3["codonpos"] = [i for i, c in enumerate(clas_entry_3["ref_codons"].split("/")[0]) if c.isupper()][0] + 1

						if str(clas_entry_3["strand"]) == "+":

							if ((
								float(clas_entry_3["codonpos"]) == 1 and 
								int(otr_entry["Start"]) >= int(clas_entry_3["Start"]) and 
								int(otr_entry["Start"]) <= int(clas_entry_3["Start"])+2
								)
								or 
								(
								float(clas_entry_3["codonpos"]) == 2 and 
								int(otr_entry["Start"]) >= int(clas_entry_3["Start"])-1 and 
								int(otr_entry["Start"]) <= int(clas_entry_3["Start"])+1
								)
								or
								(
								float(clas_entry_3["codonpos"]) == 3 and 
								int(otr_entry["Start"]) >= int(clas_entry_3["Start"])-2 and 
								int(otr_entry["Start"]) <= int(clas_entry_3["Start"])
								)):

								finalClas3SepDf = pd.concat([finalClas3SepDf, all_clas_entries_3.iloc[[z]]], ignore_index = True)
								finalOtrSepDf = pd.concat([finalOtrSepDf, otr_entry], ignore_index = True)
								
								continue

						else:

							if ((
								4-float(clas_entry_3["codonpos"]) == 1 and 
								int(otr_entry["Start"]) >= int(clas_entry_3["Start"]) and 
								int(otr_entry["Start"]) <= int(clas_entry_3["Start"])+2
								)
								or 
								(
								4-float(clas_entry_3["codonpos"]) == 2 and 
								int(otr_entry["Start"]) >= int(clas_entry_3["Start"])-1 and 
								int(otr_entry["Start"]) <= int(clas_entry_3["Start"])+1
								)
								or
								(
								4-float(clas_entry_3["codonpos"]) == 3 and 
								int(otr_entry["Start"]) >= int(clas_entry_3["Start"])-2 and 
								int(otr_entry["Start"]) <= int(clas_entry_3["Start"])
								)):

								finalClas3SepDf = pd.concat([finalClas3SepDf, all_clas_entries_3.iloc[[z]]], ignore_index = True)
								finalOtrSepDf = pd.concat([finalOtrSepDf, otr_entry], ignore_index = True)

								continue

						if str(clas_entry_3["vep_refseq"]) == otr_entry["vep_refseq"].to_string(index = False, header = False)[1:]:

							finalClas3SepDf = pd.concat([finalClas3SepDf, all_clas_entries_3.iloc[[z]]], ignore_index = True)
							finalOtrSepDf = pd.concat([finalOtrSepDf, otr_entry], ignore_index = True)

							continue

						if otr_entry["aa_sub_name"].to_string(index = False, header = False)[1:] in otrSubList:
							otrSubList.remove(otr_entry["aa_sub_name"].to_string(index = False, header = False)[1:])


			for y in indices(list(smDf["aa_sub_name"]), uniq_clas_entry["aa_sub_name"]):

				sm_entry = smDf.iloc[[y]]
				all_clas_entries_2 = clasDf[(clasDf["aa_sub_name"] == uniq_clas_entry["aa_sub_name"]) & (clasDf["simple_name"] != sm_entry["simple_name"].to_string(index = False, header = False)[1:])].reset_index(drop = True)

				for z, clas_entry_2 in all_clas_entries_2.iterrows():

					if (
						str(clas_entry_2["full_aa"]) == sm_entry["full_aa"].to_string(index = False, header = False)[1:] and 
						int(clas_entry_2["protein_position"]) == int(sm_entry["protein_position"]) and
						str(clas_entry_2["aa_sub_name"]) not in vusSubList and
						str(clas_entry_2["aa_sub_name"]) not in otrSubList
						):

						clas_entry_2["codonpos"] = [i for i, c in enumerate(clas_entry_2["ref_codons"].split("/")[0]) if c.isupper()][0] + 1

						if str(clas_entry_2["strand"]) == "+":

							if ((
								float(clas_entry_2["codonpos"]) == 1 and 
								int(sm_entry["Start"]) >= int(clas_entry_2["Start"]) and 
								int(sm_entry["Start"]) <= int(clas_entry_2["Start"])+2
								)
								or 
								(
								float(clas_entry_2["codonpos"]) == 2 and 
								int(sm_entry["Start"]) >= int(clas_entry_2["Start"])-1 and 
								int(sm_entry["Start"]) <= int(clas_entry_2["Start"])+1
								)
								or
								(
								float(clas_entry_2["codonpos"]) == 3 and 
								int(sm_entry["Start"]) >= int(clas_entry_2["Start"])-2 and 
								int(sm_entry["Start"]) <= int(clas_entry_2["Start"])
								)):

								finalClas2SepDf = pd.concat([finalClas2SepDf, all_clas_entries_2.iloc[[z]]], ignore_index = True)
								finalSmSepDf = pd.concat([finalSmSepDf, sm_entry], ignore_index = True)
								
								continue

						else:

							if ((
								4-float(clas_entry_2["codonpos"]) == 1 and 
								int(sm_entry["Start"]) >= int(clas_entry_2["Start"]) and 
								int(sm_entry["Start"]) <= int(clas_entry_2["Start"])+2
								)
								or 
								(
								4-float(clas_entry_2["codonpos"]) == 2 and 
								int(sm_entry["Start"]) >= int(clas_entry_2["Start"])-1 and 
								int(sm_entry["Start"]) <= int(clas_entry_2["Start"])+1
								)
								or
								(
								4-float(clas_entry_2["codonpos"]) == 3 and 
								int(sm_entry["Start"]) >= int(clas_entry_2["Start"])-2 and 
								int(sm_entry["Start"]) <= int(clas_entry_2["Start"])
								)):

								finalClas2SepDf = pd.concat([finalClas2SepDf, all_clas_entries_2.iloc[[z]]], ignore_index = True)
								finalSmSepDf = pd.concat([finalSmSepDf, sm_entry], ignore_index = True)

								continue

						if str(clas_entry_2["vep_refseq"]) == sm_entry["vep_refseq"].to_string(index = False, header = False)[1:]:

							finalClas2SepDf = pd.concat([finalClas2SepDf, all_clas_entries_2.iloc[[z]]], ignore_index = True)
							finalSmSepDf = pd.concat([finalSmSepDf, sm_entry], ignore_index = True)

							continue


		finalClasSepDf.drop(finalClasSepDf.filter(regex="Unname"),axis=1, inplace=True)
		finalVusSepDf.drop(finalVusSepDf.filter(regex="Unname"),axis=1, inplace=True)

		finalClas2SepDf.drop(finalClas2SepDf.filter(regex="Unname"),axis=1, inplace=True)
		finalSmSepDf.drop(finalSmSepDf.filter(regex="Unname"),axis=1, inplace=True)

		finalClas3SepDf.drop(finalClas3SepDf.filter(regex="Unname"),axis=1, inplace=True)
		finalOtrSepDf.drop(finalOtrSepDf.filter(regex="Unname"),axis=1, inplace=True)


		clas_columns_sep = []
		for column in list(finalClasSepDf.columns):
			clas_columns_sep.append(clas_v2 + "_" + str(column))
		finalClasSepDf.columns = clas_columns_sep

		vus_columns_sep = []
		for column in list(finalVusSepDf.columns):
			vus_columns_sep.append("vus_" + str(column))
		finalVusSepDf.columns = vus_columns_sep


		clas_2_columns_sep = []
		for column in list(finalClas2SepDf.columns):
			clas_2_columns_sep.append(clas_v2 + "_" + str(column))
		finalClas2SepDf.columns = clas_2_columns_sep

		sm_columns_sep = []
		for column in list(finalSmSepDf.columns):
			sm_columns_sep.append("sm_" + str(column))
		finalSmSepDf.columns = sm_columns_sep


		clas_3_columns_sep = []
		for column in list(finalClas3SepDf.columns):
			clas_3_columns_sep.append(clas_v2 + "_" + str(column))
		finalClas3SepDf.columns = clas_3_columns_sep

		otr_columns_sep = []
		for column in list(finalOtrSepDf.columns):
			otr_columns_sep.append(clas_v2_otr + "_" + str(column))
		finalOtrSepDf.columns = otr_columns_sep


		finalSepDifDf = pd.concat([finalClasSepDf, finalVusSepDf], axis = 1)
		finalSepSamDf = pd.concat([finalClas2SepDf, finalSmSepDf], axis = 1)
		finalSepOppDf = pd.concat([finalClas3SepDf, finalOtrSepDf], axis = 1)

		finalSepDifDfND = finalSepDifDf.drop_duplicates(subset = [f"{clas_v2}_aa_sub_name"]).drop_duplicates(subset = ["vus_aa_sub_name"])
		finalSepSamDfND = finalSepSamDf.drop_duplicates(subset = [f"{clas_v2}_aa_sub_name"]).drop_duplicates(subset = ["sm_aa_sub_name"])
		finalSepOppDfND = finalSepOppDf.drop_duplicates(subset = [f"{clas_v2}_aa_sub_name"]).drop_duplicates(subset = [f"{clas_v2_otr}_aa_sub_name"])
		inCommonDf = finalSepDifDf[finalSepDifDf["vus_aa_sub_name"].isin(finalSepOppDf[f"{clas_v2_otr}_aa_sub_name"])]

		finalSepDifDfND.to_csv(f"/net/data/aasubs/clinvar_only/results/strong_nd_difclas_{clas_v2}_vus.csv")
		finalSepSamDfND.to_csv(f"/net/data/aasubs/clinvar_only/results/strong_nd_difclas_{clas_v2}_{clas_v2}.csv")
		finalSepOppDfND.to_csv(f"/net/data/aasubs/clinvar_only/results/strong_nd_difclas_{clas_v2}_{clas_v2_otr}.csv")
		inCommonDf.to_csv(f"/net/data/aasubs/clinvar_only/results/strong_nd_difclas_{clas_v2}_in_common.csv")

		finalSepDifDf.sort_values(by = ["vus_aa_sub_name"]).reset_index(drop = True).to_csv(f"/net/data/aasubs/clinvar_only/results/strong_{clas_v2}_analysis.csv")

