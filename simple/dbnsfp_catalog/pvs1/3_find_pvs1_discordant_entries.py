import pandas as pd
import sys


geneList = []


def get_pvs1(cvDf, geneName):

	cvDfPLP = cvDf[cvDf["simple_annot"] == "P/LP"]

	if len(cvDf) > 0 and len(cvDfPLP)/len(cvDf) > 0.5:
		geneList.append(geneName)


def main():

	cvStopDf = pd.read_csv("/net/data/aasubs/clinvar_only/pvs1/stop_annotated_entries.csv")
	cvStopDf["id"] = "stop"
	cvSpliceDf = pd.read_csv("/net/data/aasubs/clinvar_only/pvs1/splice_annotated_entries.csv")
	cvSpliceDf["id"] = "splice"

	cvDf = pd.concat([cvStopDf, cvSpliceDf], ignore_index = True).sort_values(by = ["simple_name"]).reset_index(drop = True)

	noneDf = cvDf.fillna("-")[cvDf.fillna("-")["gene_name"] == "-"]

	cvDf = cvDf.dropna(subset = ["gene_name"])

	cvDf = cvDf.fillna("-")

	noneDf["gene_name"] = noneDf["Name"].str.split(":").str.get(0)
	noneDf["refseq"] = "-"

	cvDf = pd.concat([cvDf, noneDf], ignore_index = True).sort_values(by = ["simple_name"]).reset_index(drop = True)

	dbStopDf = pd.read_csv("/net/data/aasubs/dbnsfp_catalog/pvs1/dbnsfp_stop_annotated_entries.csv")
	dbStopDf["id"] = "stop"
	dbSpliceDf = pd.read_csv("/net/data/aasubs/dbnsfp_catalog/pvs1/dbnsfp_splice_annotated_entries.csv")
	dbSpliceDf["id"] = "splice"

	dbDf = pd.concat([dbStopDf, dbSpliceDf], ignore_index = True)

	dbDf = dbDf.fillna("-")

	dbDf = dbDf.sort_values(by = ["simple_name"]).reset_index(drop = True)

	for geneName in list(cvDf.drop_duplicates(subset = ["gene_name"])["gene_name"]):
		get_pvs1(cvDf[cvDf["gene_name"] == geneName], geneName)

	geneDf = pd.DataFrame()
	geneDf["gene_name"] = geneList

	totalDf = dbDf[dbDf["vep_gene_name"].isin(geneDf["gene_name"])]

	totalDf.drop(totalDf.filter(regex="Unname"), axis=1, inplace=True)
  totalDf = totalDf.sort_values(by = ["simple_name"]).reset_index(drop = True)

	totalDf.to_csv(f"/net/data/aasubs/dbnsfp_catalog/pvs1/dbnsfp_pvs1_discordances.csv")


if __name__ == '__main__':
	main()
