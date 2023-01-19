import pandas as pd


def main():

	pd.options.mode.chained_assignment = None

	df = pd.read_csv("/net/data/aasubs/wsplice_clinvar_11-06-22.csv")
	df.drop(df.filter(regex="Unname"), axis=1, inplace=True)

	df = df[(df["Name"].str[-5:-3] == "-1") | (df["Name"].str[-5:-3] == "-2") | (df["Name"].str[-5:-3] == "+1") | (df["Name"].str[-5:-3] == "+2")]

	df = df.sort_values(by = ["simple_name"]).reset_index(drop = True)

	df.to_csv("/net/data/aasubs/clinvar_only/pvs1/splice_entries.csv")


if __name__ == "__main__":
	main()
