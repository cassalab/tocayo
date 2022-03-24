import pandas as pd

df = pd.read_csv("/net/data/aasubs/parsed_clinvar_12-19.csv")
df.drop(df.filter(regex="Unname"), axis=1, inplace=True)

evidenceList = ["strong", "moderate"]

for evidence in evidenceList:

	aa_name = "aa_sub_name" if evidence == "strong" else "aa_half_name"

	freq = df[aa_name].value_counts()
	frequent_values = freq[freq > 1].index
	df[df[aa_name].isin(frequent_values)].sort_values(by = ["aa_sub_name"]).reset_index(drop = True).to_csv(f"/net/data/aasubs/clinvar_only/dup_files/{evidence}_all_dup_subs.csv")
