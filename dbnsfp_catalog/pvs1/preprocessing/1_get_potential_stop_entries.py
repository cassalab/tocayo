import sqlite3 
import pandas as pd
import numpy as np


def get_pvs1(cvDf, geneName):

	query = f"""SELECT chr,HGVSc_ANNOVAR,pos,genename,ref,alt,aapos,aaref,aaalt,Ensembl_transcriptid,codonpos,cds_strand FROM features WHERE (genename = "{geneName}" COLLATE NOCASE OR genename LIKE "{geneName};%")"""

	dbDf = pd.read_sql_query(query, con)

	dbDf = dbDf.dropna(subset = ["aaref", "aaalt"]).reset_index(drop = True)

	dbDf["simple_name"] = dbDf["chr"].map(str) + "-" + dbDf["pos"].map(str) + "-" + dbDf["ref"].map(str) + "-" + dbDf["alt"].map(str)

	dbDf = dbDf[(dbDf["aaalt"] == "X") & (dbDf["simple_name"].isin(cvDf["simple_name"]) == False)].reset_index(drop = True)

	dbDf["used_genename"] = geneName

	return dbDf


if __name__ == "__main__":

	global con
	con = sqlite3.connect("/net/db/dbnsfp/dbnsfp.sqlite")

	df = pd.read_csv("/net/data/aasubs/parsed_clinvar_11-06-22.csv")

	totalDf = pd.DataFrame(columns = ["chr", "HGVSc_ANNOVAR", "pos", "genename", "ref", "alt", "aapos", "aaref", "aaalt", "Ensembl_transcriptid", "codonpos", "cds_strand", "simple_name", "used_genename"])

	for geneName in list(df.drop_duplicates(subset = ["gene_name"])["gene_name"]):
		totalDf = pd.concat([totalDf, get_pvs1(df[df["gene_name"] == geneName], geneName)], ignore_index = True)

	totalDf.drop(totalDf.filter(regex="Unname"), axis=1, inplace=True)

	totalDf.sort_values(by = ["simple_name"]).reset_index(drop = True).to_csv(f"/net/data/aasubs/update_11-06-22/dbnsfp_catalog/pvs1/dbnsfp_stop_entries.csv")
