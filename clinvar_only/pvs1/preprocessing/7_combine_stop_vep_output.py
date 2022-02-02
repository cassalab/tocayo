import pandas as pd


file_starters = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", '22', "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34"]

vepDfFull = pd.DataFrame()


for x in file_starters:
	part_file_path = f"/net/data/aasubs/clinvar_only/pvs1/{x}_stop_vep_output.vcf"
	vepDfPart = pd.read_csv(part_file_path, sep = "\t", skiprows = 105)

	if x == "0":
		vepDfFull = pd.DataFrame(columns = vepDfPart.columns)

	vepDfFull = pd.concat([vepDfFull, vepDfPart], ignore_index = True)


vepDfFull.to_csv("/net/data/aasubs/clinvar_only/pvs1/stop_vep_output.csv")
