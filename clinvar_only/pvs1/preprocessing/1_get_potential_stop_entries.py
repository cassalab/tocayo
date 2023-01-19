import pandas as pd


def main():

	pd.options.mode.chained_assignment = None

	df = pd.read_csv("/net/data/aasubs/parsed_clinvar_11-06-22.csv")
	df.drop(df.filter(regex="Unname"), axis=1, inplace=True)
	df["aa_alt"] = df["aa_sub_name"].str[-3:]
	df = df[df["aa_alt"] == "Ter"]
	
	df = df.sort_values(by = ["simple_name"]).reset_index(drop = True)

	df.to_csv("/net/data/aasubs/clinvar_only/pvs1/stop_entries.csv")


if __name__ == "__main__":
	main()
