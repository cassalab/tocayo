import pandas as pd
import numpy as np
import time
import sys

pd.options.mode.chained_assignment = None

cvDf = pd.read_csv(f"/net/data/aasubs/dbnsfp_catalog/annotated/strong_annotated_cv_potential_entries.csv", dtype = str)
dbDf = pd.read_csv(f"/net/data/aasubs/dbnsfp_catalog/annotated/strong_annotated_db_potential_entries.csv", dtype = str)

cvDf.drop(cvDf.filter(regex="Unname"), axis=1, inplace=True)
dbDf.drop(dbDf.filter(regex="Unname"), axis=1, inplace=True)
dbDf["Start"] = dbDf["pos"]

cv_columns = []
for column in list(cvDf.columns):
	cv_columns.append("cv_" + column)

db_columns = []
for column in list(dbDf.columns):
	db_columns.append("db_" + column)

cvDf.columns = cv_columns
dbDf.columns = db_columns

combDf = pd.concat([cvDf, dbDf], axis = 1)
combDf = combDf.drop_duplicates(subset = ["cv_simple_name", "db_simple_name"])

outDf = combDf[((combDf["cv_impact"] == "MODERATE") | (combDf["cv_impact"] == "HIGH")) & (((combDf["db_impact"] != "MODERATE") & (combDf["db_impact"] != "HIGH")) | (combDf["cv_full_aa"] != combDf["db_full_aa"]) | (combDf["cv_protein_position"] != combDf["db_protein_position"]))]
combDf = combDf[((combDf["cv_impact"] == "MODERATE") | (combDf["cv_impact"] == "HIGH")) & ((combDf["db_impact"] == "MODERATE") | (combDf["db_impact"] == "HIGH")) & (combDf["cv_full_aa"] == combDf["db_full_aa"]) & (combDf["cv_protein_position"] == combDf["db_protein_position"])]

combDf = combDf.sort_values(by = ["cv_aa_sub_name"]).reset_index(drop = True)
outDf = outDf.sort_values(by = ["cv_aa_sub_name"]).reset_index(drop = True)

combDf_single = combDf[combDf["cv_vep_refseq"] == combDf["db_vep_refseq"]]
combDf_single.to_csv("/net/data/aasubs/dbnsfp_catalog/filtered/strong_filtered_potential_entries_single.csv")

combDf_non_single = combDf[combDf["cv_vep_refseq"] != combDf["db_vep_refseq"]].sort_values(by = ["cv_aa_sub_name"]).reset_index(drop = True)


for i, entry in combDf_non_single.iterrows():

	entry["cv_codonpos"] = [i for i, c in enumerate(entry["cv_ref_codons"].split("/")[0]) if c.isupper()][0] + 1

	if entry["cv_strand"] == "+":

		if not ((
			float(entry["cv_codonpos"]) == 1 and 
			int(entry["cv_Start"]) >= int(entry["db_Start"]) and 
			int(entry["cv_Start"]) <= int(entry["db_Start"])+2
			)
			or 
			(
			float(entry["cv_codonpos"]) == 2 and 
			int(entry["cv_Start"]) >= int(entry["db_Start"])-1 and 
			int(entry["cv_Start"]) <= int(entry["db_Start"])+1
			)
			or
			(
			float(entry["cv_codonpos"]) == 3 and 
			int(entry["cv_Start"]) >= int(entry["db_Start"])-2 and 
			int(entry["cv_Start"]) <= int(entry["db_Start"])
			)):

			combDf_non_single.drop(i, inplace = True)

	else:

		if not ((
			4-float(entry["cv_codonpos"]) == 1 and 
			int(entry["cv_Start"]) >= int(entry["db_Start"]) and 
			int(entry["cv_Start"]) <= int(entry["db_Start"])+2
			)
			or 
			(
			4-float(entry["cv_codonpos"]) == 2 and 
			int(entry["cv_Start"]) >= int(entry["db_Start"])-1 and 
			int(entry["cv_Start"]) <= int(entry["db_Start"])+1
			)
			or
			(
			4-float(entry["cv_codonpos"]) == 3 and 
			int(entry["cv_Start"]) >= int(entry["db_Start"])-2 and 
			int(entry["cv_Start"]) <= int(entry["db_Start"])
			)):

			combDf_non_single.drop(i, inplace = True)


combDf = pd.concat([combDf_single, combDf_non_single], ignore_index = True).sort_values(by = ["cv_aa_sub_name"]).reset_index(drop = True)
combDf.to_csv("/net/data/aasubs/dbnsfp_catalog/filtered/strong_filtered_potential_entries_both.csv")


vepDf = pd.read_csv("/net/data/aasubs/update_11-06-22/dbnsfp_catalog/vep/strong_db_vep_output.txt", sep = "\t", skiprows = 105, dtype = str)
vepDf["simple_name"] = vepDf["#Uploaded_variation"].str.replace("_", "-")
vepDf["simple_name"] = vepDf["simple_name"].str.replace("/", "-")
vepDf = vepDf[vepDf["simple_name"].isin(outDf["db_simple_name"])]

vepDf = vepDf.drop_duplicates(subset = ["simple_name", "Feature"])

checkDf = outDf[outDf["db_simple_name"].isin(vepDf["simple_name"])].drop_duplicates(subset = ["db_simple_name"])
vepDf = vepDf.sort_values(by = ["simple_name"]).reset_index(drop = True)
vepDf.drop(vepDf.filter(regex="Unname"), axis=1, inplace=True)

checkDfCI = outDf[outDf["db_simple_name"].isin(vepDf["simple_name"])]
checkDfCI = checkDfCI.sort_values(by = ["db_simple_name"]).reset_index(drop = True)
checkDfCI = checkDfCI.drop_duplicates(subset = ["db_simple_name"])

checkDf = checkDf.sort_values(by = ["db_simple_name"]).reset_index(drop = True)
checkDfSave = outDf[outDf["db_simple_name"].isin(vepDf["simple_name"])]
checkDf.drop(checkDf.filter(regex="Unname"), axis=1, inplace=True)

vepDfCI = vepDf.drop_duplicates(subset = ["simple_name"])
checkDf = checkDf.set_index(vepDfCI.index)
checkDf = checkDf.reindex(vepDf.index, method = "ffill").sort_values(by = ["db_simple_name"]).reset_index(drop = True)

vepFinalDf = pd.DataFrame(columns = vepDf.columns)
vepDf = vepDf[vepDf["simple_name"].isin(vepFinalDf["simple_name"]) == False]
checkDf = checkDf[checkDf["db_simple_name"].isin(vepFinalDf["simple_name"]) == False]
vepFinalDf = pd.concat([vepFinalDf, vepDf[(vepDf["CANONICAL"] == "YES") & ((vepDf["IMPACT"] == "MODERATE") | (vepDf["IMPACT"] == "HIGH")) & (vepDf["BIOTYPE"] == "protein_coding") & (vepDf["Amino_acids"] == checkDf["cv_full_aa"]) & (vepDf["Protein_position"] == checkDf["cv_protein_position"])].drop_duplicates(subset = ["simple_name"])], ignore_index = True).drop_duplicates(subset = ["simple_name"]).reset_index(drop = True)
vepDf = vepDf[vepDf["simple_name"].isin(vepFinalDf["simple_name"]) == False]
checkDf = checkDf[checkDf["db_simple_name"].isin(vepFinalDf["simple_name"]) == False]
vepFinalDf = pd.concat([vepFinalDf, vepDf[((vepDf["IMPACT"] == "MODERATE") | (vepDf["IMPACT"] == "HIGH")) & (vepDf["BIOTYPE"] == "protein_coding") & (vepDf["Amino_acids"] == checkDf["cv_full_aa"]) & (vepDf["Protein_position"] == checkDf["cv_protein_position"])].drop_duplicates(subset = ["simple_name"])], ignore_index = True).drop_duplicates(subset = ["simple_name"]).reset_index(drop = True)
vepDf = vepDf[vepDf["simple_name"].isin(vepFinalDf["simple_name"]) == False]
checkDf = checkDf[checkDf["db_simple_name"].isin(vepFinalDf["simple_name"]) == False]
vepFinalDf = pd.concat([vepFinalDf, vepDf[((vepDf["IMPACT"] == "MODERATE") | (vepDf["IMPACT"] == "HIGH")) & (vepDf["BIOTYPE"] == "protein_coding")].drop_duplicates(subset = ["simple_name"])], ignore_index = True).drop_duplicates(subset = ["simple_name"]).reset_index(drop = True)
vepDf = vepDf[vepDf["simple_name"].isin(vepFinalDf["simple_name"]) == False]
checkDf = checkDf[checkDf["db_simple_name"].isin(vepFinalDf["simple_name"]) == False]
vepFinalDf = pd.concat([vepFinalDf, vepDf[(vepDf["IMPACT"] == "LOW") & (vepDf["BIOTYPE"] == "protein_coding")].drop_duplicates(subset = ["simple_name"])], ignore_index = True).drop_duplicates(subset = ["simple_name"]).reset_index(drop = True)
vepDf = vepDf[vepDf["simple_name"].isin(vepFinalDf["simple_name"]) == False]
checkDf = checkDf[checkDf["db_simple_name"].isin(vepFinalDf["simple_name"]) == False]
vepFinalDf = pd.concat([vepFinalDf, vepDf.drop_duplicates(subset = ["simple_name"])], ignore_index = True).drop_duplicates(subset = ["simple_name"]).reset_index(drop = True)

vepFinalDf = vepFinalDf.sort_values(by = ["simple_name"]).reset_index(drop = True)

checkDfSave = checkDfSave.sort_values(by = ["db_simple_name"]).reset_index(drop = True)

vepFinalDf = vepFinalDf.set_index(checkDfCI.index)
vepFinalDf = vepFinalDf.reindex(checkDfSave.index, method = "ffill").sort_values(by = ["simple_name"]).reset_index(drop = True)


checkDfSave["db_strand"] = np.where(vepFinalDf["STRAND"] == "1", "+", "-")
checkDfSave["db_full_codons"] = vepFinalDf["Codons"].str.lower()
checkDfSave["db_half_codons"] = vepFinalDf["Codons"].str.split("/").str.get(0).str.lower()
checkDfSave["db_ref_codons"] = vepFinalDf["Codons"].str.split("/").str.get(0)
checkDfSave["db_full_aa"] = vepFinalDf["Amino_acids"]
checkDfSave["db_half_aa"] = vepFinalDf["Amino_acids"].str.split("/").str.get(0)
checkDfSave["db_protein_position"] = vepFinalDf["Protein_position"]
checkDfSave["db_consequence"] = vepFinalDf["Consequence"]
checkDfSave["db_impact"] = vepFinalDf["IMPACT"]
checkDfSave["db_1000_af"] = vepFinalDf["AF"]
checkDfSave["db_gnomade_af"] = vepFinalDf["gnomADe_AF"]
checkDfSave["db_gnomadg_af"] = vepFinalDf["gnomADg_AF"]
checkDfSave["db_max_af"] = vepFinalDf["MAX_AF"]
checkDfSave["db_cadd_raw"] = vepFinalDf["CADD_RAW"]
checkDfSave["db_cadd_phred"] = vepFinalDf["CADD_PHRED"]
checkDfSave["db_polyphen_score"] = vepFinalDf["PolyPhen"].str.split("(").str.get(1).str.split(")").str.get(0)
checkDfSave["db_sift_score"] = vepFinalDf["SIFT"].str.split("(").str.get(1).str.split(")").str.get(0)
checkDfSave["db_polyphen_clas"] = vepFinalDf["PolyPhen"].str.split("(").str.get(0)
checkDfSave["db_sift_clas"] = vepFinalDf["SIFT"].str.split("(").str.get(0)
checkDfSave["db_vep_simple_name"] = vepFinalDf["simple_name"]
checkDfSave["db_vep_gene_name"] = vepFinalDf["SYMBOL"]
checkDfSave["db_vep_refseq"] = vepFinalDf["Feature"]


checkDfSave = checkDfSave[((checkDfSave["db_impact"] == "MODERATE") | (checkDfSave["db_impact"] == "HIGH")) & (checkDfSave["db_full_aa"] == checkDfSave["cv_full_aa"]) & (checkDfSave["db_protein_position"] == checkDfSave["cv_protein_position"])]


for i, entry in checkDfSave.iterrows():

	if (
		entry["cv_vep_refseq"] != entry["db_vep_refseq"]
		):

		entry["cv_codonpos"] = [i for i, c in enumerate(entry["cv_ref_codons"].split("/")[0]) if c.isupper()][0] + 1

		if entry["cv_strand"] == "+":

			if not ((
				float(entry["cv_codonpos"]) == 1 and 
				int(entry["cv_Start"]) >= int(entry["db_Start"]) and 
				int(entry["cv_Start"]) <= int(entry["db_Start"])+2
				)
				or 
				(
				float(entry["cv_codonpos"]) == 2 and 
				int(entry["cv_Start"]) >= int(entry["db_Start"])-1 and 
				int(entry["cv_Start"]) <= int(entry["db_Start"])+1
				)
				or
				(
				float(entry["cv_codonpos"]) == 3 and 
				int(entry["cv_Start"]) >= int(entry["db_Start"])-2 and 
				int(entry["cv_Start"]) <= int(entry["db_Start"])
				)):

				checkDfSave.drop(i, inplace = True)

		else:

			if not ((
				4-float(entry["cv_codonpos"]) == 1 and 
				int(entry["cv_Start"]) >= int(entry["db_Start"]) and 
				int(entry["cv_Start"]) <= int(entry["db_Start"])+2
				)
				or 
				(
				4-float(entry["cv_codonpos"]) == 2 and 
				int(entry["cv_Start"]) >= int(entry["db_Start"])-1 and 
				int(entry["cv_Start"]) <= int(entry["db_Start"])+1
				)
				or
				(
				4-float(entry["cv_codonpos"]) == 3 and 
				int(entry["cv_Start"]) >= int(entry["db_Start"])-2 and 
				int(entry["cv_Start"]) <= int(entry["db_Start"])
				)):

				checkDfSave.drop(i, inplace = True)


combDf = pd.concat([combDf, checkDfSave], ignore_index = True).sort_values(by = ["cv_aa_sub_name"]).reset_index(drop = True)
combDf.to_csv("/net/data/aasubs/dbnsfp_catalog/filtered/strong_filtered_potential_entries_full.csv")
