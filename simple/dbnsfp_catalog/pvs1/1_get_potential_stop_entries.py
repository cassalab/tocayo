import sqlite3 
import pandas as pd
import numpy as np
import time
import sys


def get_pvs1(cvDf):

	query = f"""SELECT chr,HGVSc_ANNOVAR,pos,genename,ref,alt,aapos,aaref,aaalt,Ensembl_transcriptid,codonpos,cds_strand FROM features WHERE (genename = "{geneName}" COLLATE NOCASE OR genename LIKE "{geneName};%")"""

	dbDf = pd.read_sql_query(query, con)

	dbDf = dbDf.dropna(subset = ["aaref", "aaalt"]).reset_index(drop = True)

	dbDf["simple_name"] = dbDf["chr"].map(str) + "-" + dbDf["pos"].map(str) + "-" + dbDf["ref"].map(str) + "-" + dbDf["alt"].map(str)

	dbDf = dbDf[(dbDf["aaalt"] == "X") & (dbDf["simple_name"].isin(cvDf["simple_name"]) == False)]

	return dbDf


if __name__ == "__main__":

	global con
	con = sqlite3.connect('/net/db/dbnsfp/dbnsfp.sqlite')

	aaConv = {
		"A": "Ala",
		"R": "Arg",
		"N": "Asn",
		"D": "Asp",
		"C": "Cys",
		"E": "Glu",
		"Q": "Gln",
		"G": "Gly",
		"H": "His",
		"I": "Ile",
		"L": "Leu",
		"K": "Lys",
		"M": "Met",
		"F": "Phe",
		"P": "Pro",
		"S": "Ser",
		"T": "Thr",
		"W": "Trp",
		"Y": "Tyr",
		"V": "Val",
		"O": "Pyl",
		"U": "Sec",
		"B": "Asx",
		"Z": "Glx",
		"X": "Ter",
		"J": "Xle"
	}

	df = pd.read_csv("/net/data/aasubs/parsed_clinvar_12-19.csv")

	totalDf = pd.DataFrame(columns = ["chr", "HGVSc_ANNOVAR", "pos", "genename", "ref", "alt", "aapos", "aaref", "aaalt", "Ensembl_transcriptid", "codonpos", "cds_strand", "simple_name"])

	for geneName in list(df.drop_duplicates(subset = ["gene_name"])["gene_name"]):
		totalDf = pd.concat([totalDf, get_pvs1(df[df["gene_name"] == geneName])], ignore_index = True)

	totalDf.drop(totalDf.filter(regex="Unname"), axis=1, inplace=True)

	totalDf.sort_values(by = ["simple_name"]).reset_index(drop = True).to_csv(f"/net/data/aasubs/dbnsfp_catalog/pvs1/dbnsfp_stop_entries.csv")
