import sqlite3 
import pandas as pd
import numpy as np
import time
import sys


def get_pvs1(cvDf):

	query = f"""SELECT chr,HGVSc_ANNOVAR,pos,genename,ref,alt,aapos,aaref,aaalt,Ensembl_transcriptid,codonpos,cds_strand FROM features WHERE (genename = "{geneName}" COLLATE NOCASE OR genename LIKE "{geneName};%")"""

	dbDf = pd.read_sql_query(query, con)

	dbDf["simple_name"] = dbDf["chr"].map(str) + "-" + dbDf["pos"].map(str) + "-" + dbDf["ref"].map(str) + "-" + dbDf["alt"].map(str)

	dbDf = dbDf[(dbDf["aapos"].str.contains("-1")) & (dbDf["simple_name"].isin(cvDf["simple_name"]) == False)]

	return dbDf


if __name__ == "__main__":

	global con
	con = sqlite3.connect("/net/db/dbnsfp/dbnsfp.sqlite")

	df = pd.read_csv("/net/data/aasubs/parsed_clinvar_11-06-22.csv")

	noneDf = df.fillna("-200")[df.fillna("-200")["gene_name"] == "-200"]

	df = df.dropna(subset = ["gene_name"])

	noneDf["gene_name"] = noneDf["Name"].str.split(":").str.get(0)
	noneDf = noneDf.dropna(subset = ["gene_name"])
	noneDf["refseq"] = "not found"

	noneDf = noneDf[(noneDf["Name"].str.contains("NC_") == False) & (noneDf["Name"].str.contains("NG_") == False) & (noneDf["Name"].str.contains("NM_") == False) & (noneDf["Name"].str.contains("m") == False) & (noneDf["Name"].str.contains("c")) & (noneDf["Name"].str.contains(">"))]

	df = pd.concat([df, noneDf], ignore_index = True).sort_values(by = ["simple_name"]).reset_index(drop = True)

	totalDf = pd.DataFrame(columns = ["chr", "HGVSc_ANNOVAR", "pos", "genename", "ref", "alt", "aapos", "aaref", "aaalt", "Ensembl_transcriptid", "codonpos", "cds_strand", "simple_name"])

	for geneName in list(df.drop_duplicates(subset = ["gene_name"])["gene_name"]):
		totalDf = pd.concat([totalDf, get_pvs1(df[df["gene_name"] == geneName])], ignore_index = True)

	totalDf.drop(totalDf.filter(regex="Unname"), axis=1, inplace=True)

	totalDf.sort_values(by = ["simple_name"]).reset_index(drop = True).to_csv(f"/net/data/aasubs/dbnsfp_catalog/pvs1/dbnsfp_splice_entries.csv")

