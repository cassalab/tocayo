import sqlite3 
import pandas as pd
import numpy as np
import time
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


def get_variants(geneName):

	query = f"""SELECT chr,pos,genename,ref,alt,aapos,aaref,aaalt,Ensembl_transcriptid,codonpos,cds_strand FROM features WHERE (genename = "{geneName}" COLLATE NOCASE OR genename LIKE "{geneName};%")"""

	dbDf = pd.read_sql_query(query, con)

	dbDf = dbDf.dropna(subset = ["aaref", "aaalt"]).reset_index(drop = True)
	dbDf["simple_name"] = dbDf["chr"].map(str) + "-" + dbDf["pos"].map(str) + "-" + dbDf["ref"].map(str) + "-" + dbDf["alt"].map(str)

	cvDf_gene = cvDf[cvDf["gene_name"] == geneName]
	cvNameList = list(cvDf_gene["simple_name"])
	cvHalfList = list(cvDf_gene["aa_half_name"])
	cvSubList = list(cvDf_gene["aa_sub_name"])

	cvDfFinal = pd.DataFrame(columns = cvDf_gene.columns)
	dbDfFinal = pd.DataFrame(columns = dbColumns)

	for i, db_entry in dbDf.iterrows():
		for j, _ in enumerate(db_entry["Ensembl_transcriptid"].split(";")):
			simple_name = str(db_entry["simple_name"])
			aa_half_name = db_entry["genename"].split(";")[j] + "-" + aaConv.get(db_entry["aaref"]) + db_entry["aapos"].split(";")[j]
			aa_sub_name = str(db_entry["genename"]) + "-" + aaConv.get(db_entry["aaref"]) + str(db_entry["aapos"].split(";")[j]) + aaConv.get(db_entry["aaalt"])
			if aa_half_name in cvHalfList and simple_name not in cvNameList:
				for y in indices(cvHalfList, aa_half_name):
					if aa_sub_name != cvSubList[y]:
						cvDfFinal = pd.concat([cvDfFinal, cvDf_gene.iloc[[y]]], ignore_index = True)
						dbDfFinal = pd.concat([dbDfFinal, dbDf.iloc[[i]]], ignore_index = True)
				break

	return [cvDfFinal, dbDfFinal]


if __name__ == "__main__":

	global con
	con = sqlite3.connect("/net/db/dbnsfp/dbnsfp.sqlite")

	cvDf = pd.read_csv("/net/data/aasubs/parsed_clinvar_11-06-22.csv")

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

	geneList = list(cvDf.drop_duplicates(subset = ["gene_name"])["gene_name"])

	dbColumns = ["chr", "pos", "genename", "ref", "alt", "aapos", "aaref", "aaalt", "Ensembl_transcriptid", "codonpos", "cds_strand", "simple_name"]

	cvDfFull = pd.DataFrame(columns = cvDf.columns)
	dbDfFull = pd.DataFrame(columns = dbColumns)

	for geneName in geneList:
		dfs = get_variants(geneName)
		cvDfFull = pd.concat([cvDfFull, dfs[0]], ignore_index = True)
		dbDfFull = pd.concat([dbDfFull, dfs[1]], ignore_index = True)

	cvDfFull.to_csv("/net/data/aasubs/dbnsfp_catalog/moderate_cv_potential_entries.csv")
	dbDfFull.to_csv("/net/data/aasubs/dbnsfp_catalog/moderate_db_potential_entries.csv")
