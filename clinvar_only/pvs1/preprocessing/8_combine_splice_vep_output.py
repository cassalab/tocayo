import pandas as pd


file_starters = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18"]

vepDfFull = pd.DataFrame()


for x in file_starters:
	part_file_path = f"/net/data/aasubs/clinvar_only/pvs1/{x}_splice_vep_output.vcf"
	vepDfPart = pd.read_csv(part_file_path, sep = "\t", skiprows = 105)

	if x == "0":
		vepDfFull = pd.DataFrame(columns = vepDfPart.columns)

	vepDfFull = pd.concat([vepDfFull, vepDfPart], ignore_index = True)


vepDfFull.to_csv("/net/data/aasubs/clinvar_only/pvs1/splice_vep_output.csv")
