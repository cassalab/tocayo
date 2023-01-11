import pandas as pd
import numpy as np
import time
import sys


pd.options.mode.chained_assignment = None

evidenceList = ["strong", "moderate"]


def main():

	for evidence in evidenceList:

		if evidence == "strong":
			vepDf = pd.read_csv("/net/data/aasubs/dbnsfp_catalog/vep/strong_db_vep_output.txt", sep = "\t", skiprows = 105, dtype = str)
		else:
			vepDf = pd.read_csv("/net/data/aasubs/dbnsfp_catalog/vep/moderate_db_vep_output.csv", dtype = str)

		vepDf["simple_name"] = vepDf["#Uploaded_variation"].str.replace("_", "-")
		vepDf["simple_name"] = vepDf["simple_name"].str.replace("/", "-")
		vepDf = vepDf.sort_values(by = ["simple_name"]).reset_index(drop = True)
		vepDf.drop(vepDf.filter(regex="Unname"), axis=1, inplace=True)

		dbDf = pd.read_csv(f"/net/data/aasubs/dbnsfp_catalog/{evidence}_db_potential_entries.csv")
		cvDf = pd.read_csv(f"/net/data/aasubs/dbnsfp_catalog/{evidence}_cv_potential_entries.csv")
		dbDf["refseq"] = cvDf["refseq"]
		dbDf["aa_sub_name"] = cvDf["aa_sub_name"]

		dbDf.drop(dbDf.filter(regex="Unname"), axis=1, inplace=True)
		dbDf = dbDf.sort_values(by = ["simple_name"]).reset_index(drop = True)

		dbDfND = dbDf
		dbDfCI = dbDf.drop_duplicates(subset = ["simple_name"])
		dbDf = dbDf.drop_duplicates(subset = ["simple_name"]).reset_index(drop = True)

		dbDfVep = pd.DataFrame()

		for column in list(vepDf.columns):
			if column == "Feature":
				dbDfVep[column] = dbDf["refseq"]
			elif column == "simple_name":
				dbDfVep[column] = dbDf[column]
			else:
				dbDfVep[column] = "-"

		addVepDf = vepDf[vepDf["Feature"].str.contains("NM_")]
		addVepDf = addVepDf.drop_duplicates(subset = ["simple_name", "Feature"])

		combDf = pd.concat([addVepDf, dbDfVep], ignore_index = True)

		combFilDf = combDf.groupby(["Feature", "simple_name"], as_index = False).first()
		sizeDf = combDf.groupby(["Feature", "simple_name"], as_index = False).size()
		combFilDf["size"] = sizeDf["size"]
		combFilDf = combFilDf[combFilDf["size"] > 1].drop_duplicates(subset = ["simple_name"])
		combFilDf = combFilDf.sort_values(by = ["simple_name"]).reset_index(drop = True)

		vepFinalDf = combFilDf
		vepDf = vepDf[vepDf["simple_name"].isin(vepFinalDf["simple_name"]) == False]
		vepFinalDf = pd.concat([vepFinalDf, vepDf[(vepDf["CANONICAL"] == "YES") & ((vepDf["IMPACT"] == "MODERATE") | (vepDf["IMPACT"] == "HIGH"))].drop_duplicates(subset = ["simple_name"])], ignore_index = True).drop_duplicates(subset = ["simple_name"]).reset_index(drop = True)
		vepDf = vepDf[vepDf["simple_name"].isin(vepFinalDf["simple_name"]) == False]
		vepFinalDf = pd.concat([vepFinalDf, vepDf[((vepDf["IMPACT"] == "MODERATE") | (vepDf["IMPACT"] == "HIGH"))].drop_duplicates(subset = ["simple_name"])], ignore_index = True).drop_duplicates(subset = ["simple_name"]).reset_index(drop = True)
		vepDf = vepDf[vepDf["simple_name"].isin(vepFinalDf["simple_name"]) == False]
		vepFinalDf = pd.concat([vepFinalDf, vepDf[(vepDf["IMPACT"] == "LOW")].drop_duplicates(subset = ["simple_name"])], ignore_index = True).drop_duplicates(subset = ["simple_name"]).reset_index(drop = True)
		vepDf = vepDf[vepDf["simple_name"].isin(vepFinalDf["simple_name"]) == False]
		vepFinalDf = pd.concat([vepFinalDf, vepDf.drop_duplicates(subset = ["simple_name"])], ignore_index = True).drop_duplicates(subset = ["simple_name"]).reset_index(drop = True)

		vepFinalDf = vepFinalDf.sort_values(by = ["simple_name"]).reset_index(drop = True)

		dbDf["strand"] = np.where(vepFinalDf["STRAND"] == "1", "+", "-")
		dbDf["full_codons"] = vepFinalDf["Codons"].str.lower()
		dbDf["half_codons"] = vepFinalDf["Codons"].str.split("/").str.get(0).str.lower()
		dbDf["ref_codons"] = vepFinalDf["Codons"].str.split("/").str.get(0)
		dbDf["full_aa"] = vepFinalDf["Amino_acids"]
		dbDf["half_aa"] = vepFinalDf["Amino_acids"].str.split("/").str.get(0)
		dbDf["protein_position"] = vepFinalDf["Protein_position"]
		dbDf["consequence"] = vepFinalDf["Consequence"]
		dbDf["impact"] = vepFinalDf["IMPACT"]
		dbDf["1000_af"] = vepFinalDf["AF"]
		dbDf["gnomade_af"] = vepFinalDf["gnomADe_AF"]
		dbDf["gnomadg_af"] = vepFinalDf["gnomADg_AF"]
		dbDf["max_af"] = vepFinalDf["MAX_AF"]
		dbDf["cadd_raw"] = vepFinalDf["CADD_RAW"]
		dbDf["cadd_phred"] = vepFinalDf["CADD_PHRED"]
		dbDf["polyphen_score"] = vepFinalDf["PolyPhen"].str.split("(").str.get(1).str.split(")").str.get(0)
		dbDf["sift_score"] = vepFinalDf["SIFT"].str.split("(").str.get(1).str.split(")").str.get(0)
		dbDf["polyphen_clas"] = vepFinalDf["PolyPhen"].str.split("(").str.get(0)
		dbDf["sift_clas"] = vepFinalDf["SIFT"].str.split("(").str.get(0)
		dbDf["vep_simple_name"] = vepFinalDf["simple_name"]
		dbDf["vep_gene_name"] = vepFinalDf["SYMBOL"]
		dbDf["vep_refseq"] = vepFinalDf["Feature"]

		dbDf = dbDf.set_index(dbDfCI.index)
		dbDf = dbDf.reindex(dbDfND.index, method = "ffill").sort_values(by = ["aa_sub_name"]).reset_index(drop = True)

		dbDf.to_csv(f"/net/data/aasubs/dbnsfp_catalog/annotated/{evidence}_annotated_db_potential_entries.csv")


if __name__ == "__main__":
	main()
