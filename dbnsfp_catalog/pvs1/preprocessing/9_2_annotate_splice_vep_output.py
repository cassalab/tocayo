import pandas as pd
import numpy as np
import time
import sys


pd.options.mode.chained_assignment = None


def main():

	vepDf = pd.read_csv("/net/data/aasubs/dbnsfp_catalog/pvs1/dbnsfp_splice_vep_output.csv", dtype = str)

	vepDf["simple_name"] = vepDf["#Uploaded_variation"].str.replace("_", "-")
	vepDf["simple_name"] = vepDf["simple_name"].str.replace("/", "-")
	vepDf = vepDf.sort_values(by = ["simple_name"]).reset_index(drop = True)
	vepDf.drop(vepDf.filter(regex="Unname"), axis=1, inplace=True)

	dbDf = pd.read_csv(f"/net/data/aasubs/dbnsfp_catalog/pvs1/dbnsfp_splice_entries.csv")
	dbDf = dbDf.sort_values(by = ["simple_name"]).reset_index(drop = True)
	dbDf.drop(dbDf.filter(regex="Unname"), axis=1, inplace=True)

	vepFinalDf = vepDf[(vepDf["CANONICAL"] == "YES") & ((vepDf["Consequence"].str.contains("splice_acceptor_variant")) | (vepDf["Consequence"].str.contains("splice_donor_variant"))) & (vepDf["BIOTYPE"] == "protein_coding")].drop_duplicates(subset = ["simple_name"]).reset_index(drop = True)
	vepDf = vepDf[vepDf["simple_name"].isin(vepFinalDf["simple_name"]) == False]
	vepFinalDf = pd.concat([vepFinalDf, vepDf[((vepDf["Consequence"].str.contains("splice_acceptor_variant")) | (vepDf["Consequence"].str.contains("splice_donor_variant"))) & (vepDf["BIOTYPE"] == "protein_coding")].drop_duplicates(subset = ["simple_name"])], ignore_index = True)
	vepDf = vepDf[vepDf["simple_name"].isin(vepFinalDf["simple_name"]) == False]
	vepFinalDf = pd.concat([vepFinalDf, vepDf[(vepDf["Consequence"].str.contains("splice_acceptor_variant")) | (vepDf["Consequence"].str.contains("splice_donor_variant"))].drop_duplicates(subset = ["simple_name"])], ignore_index = True)
	vepDf = vepDf[vepDf["simple_name"].isin(vepFinalDf["simple_name"]) == False]
	vepFinalDf = pd.concat([vepFinalDf, vepDf.drop_duplicates(subset = ["simple_name"])], ignore_index = True)

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
	dbDf["gnomad_af"] = vepFinalDf["gnomAD_AF"]
	dbDf["max_af"] = vepFinalDf["MAX_AF"]
	dbDf["cadd_raw"] = vepFinalDf["CADD_RAW"]
	dbDf["cadd_phred"] = vepFinalDf["CADD_PHRED"]
	dbDf["vep_simple_name"] = vepFinalDf["simple_name"]
	dbDf["vep_gene_name"] = vepFinalDf["SYMBOL"]
	dbDf["vep_refseq"] = vepFinalDf["Feature"]
	dbDf["given_ref"] = vepFinalDf["GIVEN_REF"]
	dbDf["used_ref"] = vepFinalDf["USED_REF"]
	dbDf["lof_clas"] = vepFinalDf["LoF"]
	dbDf["lof_filter"] = vepFinalDf["LoF_filter"]
	dbDf["lof_flags"] = vepFinalDf["LoF_flags"]
	dbDf["lof_info"] = vepFinalDf["LoF_info"]
	dbDf["biotype"] = vepFinalDf["BIOTYPE"]

	dbDf = dbDf[(dbDf["consequence"].str.contains("splice_acceptor_variant")) | (dbDf["consequence"].str.contains("splice_donor_variant"))].sort_values(by = ["simple_name"]).reset_index(drop = True)

	dbDf.to_csv(f"/net/data/aasubs/dbnsfp_catalog/pvs1/dbnsfp_splice_annotated_entries.csv")


if __name__ == "__main__":
	main()

