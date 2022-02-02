import pandas as pd


def get_pvs1(cvDf):

	cvDf["aa_alt"] = cvDf["aa_sub_name"].str[-3:]
	cvDf = cvDf[cvDf["aa_alt"] == "Ter"]
	cvDfPLP = cvDf[cvDf["simple_annot"] == "P/LP"]

	if len(cvDf) > 0 and len(cvDfPLP)/len(cvDf) > 0.5:
		return cvDf
	else:
		return pd.DataFrame()


def main():

	pd.options.mode.chained_assignment = None

	df = pd.read_csv("/net/data/aasubs/parsed_clinvar_12-19.csv")
	df.drop(df.filter(regex="Unname"), axis=1, inplace=True)
	df["aa_alt"] = df["aa_sub_name"].str[-3:]
	df = df[df["aa_alt"] == "Ter"]

	df.to_csv("/net/data/aasubs/clinvar_only/pvs1/stop_entries.csv")


if __name__ == '__main__':
	main()

