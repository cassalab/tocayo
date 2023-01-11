import pandas as pd
import numpy as np
import time
import sys


pd.options.mode.chained_assignment = None

evidenceList = ["strong", "moderate"]


def main():

	for evidence in evidenceList:

		if evidence == "strong":
			vepDf = pd.read_csv("/net/data/aasubs/dbnsfp_catalog/vep/strong_cv_vep_output.txt", sep = "\t", skiprows = 105, dtype = str)
		else:
			vepDf = pd.read_csv("/net/data/aasubs/dbnsfp_catalog/vep/moderate_cv_vep_output.csv", dtype = str)

		vepDf["simple_name"] = vepDf["#Uploaded_variation"].str.replace("_", "-")
		vepDf["simple_name"] = vepDf["simple_name"].str.replace("/", "-")
		vepDf = vepDf.sort_values(by = ["simple_name"]).reset_index(drop = True)
		vepDf.drop(vepDf.filter(regex="Unname"), axis=1, inplace=True)

		cvDf = pd.read_csv(f"/net/data/aasubs/dbnsfp_catalog/{evidence}_cv_potential_entries.csv").sort_values(by = ["simple_name"]).reset_index(drop = True)
		cvDf.drop(cvDf.filter(regex="Unname"), axis=1, inplace=True)

		cvDfND = cvDf
		cvDfCI = cvDf.drop_duplicates(subset = ["simple_name"])
		cvDf = cvDf.drop_duplicates(subset = ["simple_name"]).reset_index(drop = True)

		cvDfVep = pd.DataFrame()

		for column in list(vepDf.columns):
			if column == "Feature":
				cvDfVep[column] = cvDf["refseq"]
			elif column == "simple_name":
				cvDfVep[column] = cvDf[column]
			else:
				cvDfVep[column] = "-"

		addVepDf = vepDf[vepDf["Feature"].str.contains("NM_")]
		addVepDf = addVepDf.drop_duplicates(subset = ["simple_name", "Feature"])

		combDf = pd.concat([addVepDf, cvDfVep], ignore_index = True)

		combFilDf = combDf.groupby(["Feature", "simple_name"], as_index = False).first()
		sizeDf = combDf.groupby(["Feature", "simple_name"], as_index = False).size()
		combFilDf["size"] = sizeDf["size"]
		combFilDf = combFilDf[combFilDf["size"] > 1].drop_duplicates(subset = ["simple_name"])
		combFilDf = combFilDf.sort_values(by = ["simple_name"]).reset_index(drop = True)
		
		vepFinalDf = combFilDf
		vepDf = vepDf[vepDf["simple_name"].isin(vepFinalDf["simple_name"]) == False]
		vepFinalDf = pd.concat([vepFinalDf, vepDf[(vepDf["CANONICAL"] == "YES") & ((vepDf["IMPACT"] == "MODERATE") | (vepDf["IMPACT"] == "HIGH")) & (vepDf["BIOTYPE"] == "protein_coding")].drop_duplicates(subset = ["simple_name"])], ignore_index = True).drop_duplicates(subset = ["simple_name"]).reset_index(drop = True)
		vepDf = vepDf[vepDf["simple_name"].isin(vepFinalDf["simple_name"]) == False]
		vepFinalDf = pd.concat([vepFinalDf, vepDf[((vepDf["IMPACT"] == "MODERATE") | (vepDf["IMPACT"] == "HIGH")) & (vepDf["BIOTYPE"] == "protein_coding")].drop_duplicates(subset = ["simple_name"])], ignore_index = True).drop_duplicates(subset = ["simple_name"]).reset_index(drop = True)
		vepDf = vepDf[vepDf["simple_name"].isin(vepFinalDf["simple_name"]) == False]
		vepFinalDf = pd.concat([vepFinalDf, vepDf[(vepDf["IMPACT"] == "LOW") & (vepDf["BIOTYPE"] == "protein_coding")].drop_duplicates(subset = ["simple_name"])], ignore_index = True).drop_duplicates(subset = ["simple_name"]).reset_index(drop = True)
		vepDf = vepDf[vepDf["simple_name"].isin(vepFinalDf["simple_name"]) == False]
		vepFinalDf = pd.concat([vepFinalDf, vepDf.drop_duplicates(subset = ["simple_name"])], ignore_index = True).drop_duplicates(subset = ["simple_name"]).reset_index(drop = True)

		vepFinalDf = vepFinalDf.sort_values(by = ["simple_name"]).reset_index(drop = True)

		cvDf["strand"] = np.where(vepFinalDf["STRAND"] == "1", "+", "-")
		cvDf["full_codons"] = vepFinalDf["Codons"].str.lower()
		cvDf["half_codons"] = vepFinalDf["Codons"].str.split("/").str.get(0).str.lower()
		cvDf["ref_codons"] = vepFinalDf["Codons"].str.split("/").str.get(0)
		cvDf["full_aa"] = vepFinalDf["Amino_acids"]
		cvDf["half_aa"] = vepFinalDf["Amino_acids"].str.split("/").str.get(0)
		cvDf["protein_position"] = vepFinalDf["Protein_position"]
		cvDf["consequence"] = vepFinalDf["Consequence"]
		cvDf["impact"] = vepFinalDf["IMPACT"]
		cvDf["1000_af"] = vepFinalDf["AF"]
		cvDf["gnomad_af"] = vepFinalDf["gnomAD_AF"]
		cvDf["max_af"] = vepFinalDf["MAX_AF"]
		cvDf["cadd_raw"] = vepFinalDf["CADD_RAW"]
		cvDf["cadd_phred"] = vepFinalDf["CADD_PHRED"]
		cvDf["polyphen_score"] = vepFinalDf["PolyPhen"].str.split("(").str.get(1).str.split(")").str.get(0)
		cvDf["sift_score"] = vepFinalDf["SIFT"].str.split("(").str.get(1).str.split(")").str.get(0)
		cvDf["polyphen_clas"] = vepFinalDf["PolyPhen"].str.split("(").str.get(0)
		cvDf["sift_clas"] = vepFinalDf["SIFT"].str.split("(").str.get(0)
		cvDf["vep_simple_name"] = vepFinalDf["simple_name"]
		cvDf["vep_gene_name"] = vepFinalDf["SYMBOL"]
		cvDf["vep_refseq"] = vepFinalDf["Feature"]

		cvDf = cvDf.set_index(cvDfCI.index)
		cvDf = cvDf.reindex(cvDfND.index, method = "ffill").sort_values(by = ["aa_sub_name"]).reset_index(drop = True)

		cvDf.to_csv(f"/net/data/aasubs/dbnsfp_catalog/annotated/{evidence}_annotated_cv_potential_entries.csv")


if __name__ == "__main__":
	main()
