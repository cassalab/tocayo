import pandas as pd
import sys


def get_pvs1(cvDf):

	cvDfHC = cvDf[cvDf["lof_clas"] == "HC"]

	cvDfPLP = cvDfHC[cvDfHC["simple_annot"] == "P/LP"]

	if len(cvDfHC) > 0 and len(cvDfPLP)/len(cvDfHC) > 0.5:
		return cvDf
	else:
		return pd.DataFrame()


def main():

	stopDf = pd.read_csv("/net/data/aasubs/clinvar_only/pvs1/stop_annotated_entries.csv")
	stopDf["id"] = "stop"
	spliceDf = pd.read_csv("/net/data/aasubs/clinvar_only/pvs1/splice_annotated_entries.csv")
	spliceDf["id"] = "splice"

	df = pd.concat([stopDf, spliceDf], ignore_index = True).sort_values(by = ["simple_name"]).reset_index(drop = True)

	for geneName in list(df.drop_duplicates(subset = ["gene_name"])["gene_name"]):
		totalDf = pd.concat([totalDf, get_pvs1(df[df["gene_name"] == geneName])], ignore_index = True)

	totalDf.drop(totalDf.filter(regex="Unname"), axis=1, inplace=True)
	totalDf = totalDf[totalDf["simple_annot"] == "VUS"].sort_values(by = ["simple_name"]).reset_index(drop = True)

	totalDf.to_csv(f"/net/data/aasubs/clinvar_only/pvs1/pvs1_discordances.csv")


if __name__ == "__main__":
	main()
