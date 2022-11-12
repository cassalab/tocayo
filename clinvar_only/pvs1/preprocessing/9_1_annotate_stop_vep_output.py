import pandas as pd
import numpy as np
import time
import sys


pd.options.mode.chained_assignment = None


def main():

	vepDf = pd.read_csv("/net/data/aasubs/clinvar_only/pvs1/stop_vep_output.csv", dtype = str)

	vepDf["simple_name"] = vepDf["#Uploaded_variation"].str.replace("_", "-")
	vepDf["simple_name"] = vepDf["simple_name"].str.replace("/", "-")
	vepDf = vepDf.sort_values(by = ["simple_name"]).reset_index(drop = True)
	vepDf.drop(vepDf.filter(regex="Unname"), axis=1, inplace=True)

	cvDf = pd.read_csv(f"/net/data/aasubs/clinvar_only/pvs1/stop_entries.csv").sort_values(by = ["simple_name"]).reset_index(drop = True)
	cvDf.drop(cvDf.filter(regex="Unname"), axis=1, inplace=True)

	cvDfVep = pd.DataFrame()

	for column in list(vepDf.columns):
		if column == "Feature":
			cvDfVep[column] = cvDf["refseq"]
		elif column == "simple_name":
			cvDfVep[column] = cvDf[column]
		else:
			cvDfVep[column] = "-"

	combDf = pd.concat([vepDf, cvDfVep], ignore_index = True)

	combFilDf = combDf.groupby(["Feature", "simple_name"], as_index = False).first()
	sizeDf = combDf.groupby(["Feature", "simple_name"], as_index = False).size()
	combFilDf["size"] = sizeDf["size"]
	combFilDf = combFilDf[combFilDf["size"] > 1]
	combFilDf = combFilDf.sort_values(by = ["simple_name"]).reset_index(drop = True)

	vepFinalDf = combFilDf
	vepDf = vepDf[vepDf["simple_name"].isin(vepFinalDf["simple_name"]) == False]
	vepFinalDf = pd.concat([vepFinalDf, vepDf[(vepDf["CANONICAL"] == "YES") & (vepDf["Consequence"].str.contains("stop_gained"))].drop_duplicates(subset = ["simple_name"])], ignore_index = True)
	vepDf = vepDf[vepDf["simple_name"].isin(vepFinalDf["simple_name"]) == False]
	vepFinalDf = pd.concat([vepFinalDf, vepDf[vepDf["Consequence"].str.contains("stop_gained")].drop_duplicates(subset = ["simple_name"])], ignore_index = True)
	vepDf = vepDf[vepDf["simple_name"].isin(vepFinalDf["simple_name"]) == False]
	vepFinalDf = pd.concat([vepFinalDf, vepDf.drop_duplicates(subset = ["simple_name"])], ignore_index = True)

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
	cvDf["gnomad_af"] = vepFinalDf["gnomADe_AF"]
	cvDf["max_af"] = vepFinalDf["MAX_AF"]
	cvDf["cadd_raw"] = vepFinalDf["CADD_RAW"]
	cvDf["cadd_phred"] = vepFinalDf["CADD_PHRED"]
	cvDf["vep_simple_name"] = vepFinalDf["simple_name"]
	cvDf["vep_gene_name"] = vepFinalDf["SYMBOL"]
	cvDf["vep_refseq"] = vepFinalDf["Feature"]
	cvDf["lof_clas"] = vepFinalDf["LoF"]
	cvDf["lof_filter"] = vepFinalDf["LoF_filter"]
	cvDf["lof_flags"] = vepFinalDf["LoF_flags"]
	cvDf["lof_info"] = vepFinalDf["LoF_info"]

	cvDf = cvDf[cvDf["consequence"].str.contains("stop_gained")].sort_values(by = ["simple_name"]).reset_index(drop = True)

	cvDf.to_csv(f"/net/data/aasubs/clinvar_only/pvs1/stop_annotated_entries.csv")


if __name__ == '__main__':
	main()

