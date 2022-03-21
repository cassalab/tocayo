import pandas as pd
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


if __name__ == "__main__":

	pd.options.mode.chained_assignment = None

	df = pd.read_csv(f"/net/data/aasubs/clinvar_only/annotated/strong_annotated_dup_subs.csv")
	df = df[(df["impact"] == "MODERATE") | (df["impact"] == "HIGH")]
	df.drop(df.filter(regex="Unname"),axis=1, inplace=True)
	df = df.reset_index(drop = True)

	ndDf = df.drop_duplicates(subset = ["aa_sub_name"]).reset_index(drop = True)

	nameList = []

	for _, uniq_entry in ndDf.iterrows():

		for y in indices(list(df["aa_sub_name"]), uniq_entry["aa_sub_name"]):

			entry = df.iloc[[y]]
			all_entries = df[df["aa_sub_name"] == uniq_entry["aa_sub_name"]].reset_index(drop = True)

			for z, entry_2 in all_entries.iterrows():

				if (
					str(entry_2["full_aa"]) == entry["full_aa"].to_string(index = False, header = False)[1:] and 
					int(entry_2["protein_position"]) == int(entry["protein_position"]) and
					str(entry_2["simple_name"]) != df.iloc[[y]]["simple_name"].to_string(index = False, header = False)[1:]
					):

					entry_2["codonpos"] = [i for i, c in enumerate(entry_2["ref_codons"].split("/")[0]) if c.isupper()][0] + 1

					if str(entry_2["strand"]) == "+":

						if ((
							float(entry_2["codonpos"]) == 1 and 
							int(entry["Start"]) >= int(entry_2["Start"]) and 
							int(entry["Start"]) <= int(entry_2["Start"])+2
							)
							or 
							(
							float(entry_2["codonpos"]) == 2 and 
							int(entry["Start"]) >= int(entry_2["Start"])-1 and 
							int(entry["Start"]) <= int(entry_2["Start"])+1
							)
							or
							(
							float(entry_2["codonpos"]) == 3 and 
							int(entry["Start"]) >= int(entry_2["Start"])-2 and 
							int(entry["Start"]) <= int(entry_2["Start"])
							)):
							
							nameList.append(entry_2["simple_name"])
							nameList.append(entry["simple_name"].to_string(index = False, header = False)[1:])
							
							continue

					else:

						if ((
							4-float(entry_2["codonpos"]) == 1 and 
							int(entry["Start"]) >= int(entry_2["Start"]) and 
							int(entry["Start"]) <= int(entry_2["Start"])+2
							)
							or 
							(
							4-float(entry_2["codonpos"]) == 2 and 
							int(entry["Start"]) >= int(entry_2["Start"])-1 and 
							int(entry["Start"]) <= int(entry_2["Start"])+1
							)
							or
							(
							4-float(entry_2["codonpos"]) == 3 and 
							int(entry["Start"]) >= int(entry_2["Start"])-2 and 
							int(entry["Start"]) <= int(entry_2["Start"])
							)):

							nameList.append(entry_2["simple_name"])
							nameList.append(entry["simple_name"].to_string(index = False, header = False)[1:])

							continue

					if str(entry_2["vep_refseq"]) == entry["vep_refseq"].to_string(index = False, header = False)[1:]:

						nameList.append(entry_2["simple_name"])
						nameList.append(entry["simple_name"].to_string(index = False, header = False)[1:])

						continue

	nameDf = pd.DataFrame()
	nameDf["name"] = nameList
	nameDf = nameDf.drop_duplicates(subset = ["name"]).sort_values(by = ["name"]).reset_index(drop = True)
	nameDf.to_csv("/net/data/aasubs/clinvar_only/dup_files/ps1_all_subs.csv")

