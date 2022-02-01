import pandas as pd


file_starters = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "22", "23", "24"]

cvDf = pd.read_csv("/net/data/aasubs/clinvar_only/dup_files/moderate_all_dup_subs.csv")

vepDfFull = pd.DataFrame()


for x in file_starters:
	part_file_path = f"/net/data/aasubs/clinvar_only/vep/{x}_moderate_vep_output.txt"
	vepDfPart = pd.read_csv(part_file_path, sep = "\t", skiprows = 102)

	if x == "0":
		vepDfFull = pd.DataFrame(columns = vepDfPart.columns)

	vepDfFull = pd.concat([vepDfFull, vepDfPart], ignore_index = True)


vepDfFull.to_csv("/net/data/aasubs/clinvar_only/dup_files/moderate_vep_output.csv")
