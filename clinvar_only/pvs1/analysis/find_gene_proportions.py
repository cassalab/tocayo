import pandas as pd
import sys


proportion_list = []


def get_pvs1(cvDf):

	cvDfPLP = cvDf[cvDf["simple_annot"] == "P/LP"]

	proportion_list.append(len(cvDfPLP)/len(cvDf))


def main():

	stopDf = pd.read_csv("/net/data/aasubs/clinvar_only/pvs1/stop_annotated_entries.csv")
	stopDf["id"] = "stop"
	spliceDf = pd.read_csv("/net/data/aasubs/clinvar_only/pvs1/splice_annotated_entries.csv")
	spliceDf["id"] = "splice"

	df = pd.concat([stopDf, spliceDf], ignore_index = True).sort_values(by = ["simple_name"]).reset_index(drop = True)

	df = df[df["lof_clas"] == "HC"].sort_values(by = ["simple_name"]).reset_index(drop = True)

	noneDf = df.fillna("-")[df.fillna("-")["gene_name"] == "-"]

	df = df.dropna(subset = ["gene_name"])

	df = df.fillna("-")

	noneDf["gene_name"] = noneDf["Name"].str.split(":").str.get(0)
	noneDf["refseq"] = "-"

	df = pd.concat([df, noneDf], ignore_index = True).sort_values(by = ["simple_name"]).reset_index(drop = True)

	for geneName in list(df.drop_duplicates(subset = ["gene_name"])["gene_name"]):
		get_pvs1(df[df["gene_name"] == geneName])

	totalDf = pd.DataFrame()
	totalDf["proportion"] = proportion_list

	totalDf.to_csv(f"/net/data/aasubs/clinvar_only/pvs1/pvs1_proportions.csv")


if __name__ == '__main__':
	main()

