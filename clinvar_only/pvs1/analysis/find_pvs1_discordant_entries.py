import pandas as pd
import sys


def get_pvs1(cvDf):

	cvDfPLP = cvDf[cvDf["simple_annot"] == "P/LP"]

	if len(cvDfPLP)/len(cvDf) > 0.8:
		return cvDf
	else:
		return pd.DataFrame()


def main():

	stopDf = pd.read_csv("/net/data/aasubs/clinvar_only/pvs1/stop_annotated_entries.csv")
	stopDf["id"] = "stop"
	spliceDf = pd.read_csv("/net/data/aasubs/clinvar_only/pvs1/splice_annotated_entries.csv")
	spliceDf["id"] = "splice"

	df = pd.concat([stopDf, spliceDf], ignore_index = True).sort_values(by = ["simple_name"]).reset_index(drop = True)

	totalDf = pd.DataFrame(columns = df.columns)

	noneDf = df.fillna("-200")[df.fillna("-200")["gene_name"] == "-200"]

	df = df.dropna(subset = ["gene_name"])

	noneDf["gene_name"] = noneDf["Name"].str.split(":").str.get(0)

	df = pd.concat([df, noneDf], ignore_index = True).sort_values(by = ["simple_name"]).reset_index(drop = True)

	for geneName in list(df.drop_duplicates(subset = ["gene_name"])["gene_name"]):
		totalDf = pd.concat([totalDf, get_pvs1(df[df["gene_name"] == geneName])], ignore_index = True)

	totalDf.drop(totalDf.filter(regex="Unname"), axis=1, inplace=True)

	totalDf.sort_values(by = ["aa_sub_name"]).reset_index(drop = True).to_csv(f"/net/data/aasubs/clinvar_only/pvs1/pvs1_all_matches_08.csv")


if __name__ == '__main__':
	main()
