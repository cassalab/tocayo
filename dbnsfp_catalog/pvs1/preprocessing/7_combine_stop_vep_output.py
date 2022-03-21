import pandas as pd


file_starters = []

for i in range(124):
	file_starters.append(f"{i}")

vepDfFull = pd.DataFrame()


for x in file_starters:
	part_file_path = f"/net/data/aasubs/dbnsfp_catalog/pvs1/{x}_dbn_stop_vep_output.vcf"
	vepDfPart = pd.read_csv(part_file_path, sep = "\t", skiprows = 105)

	if x == "0":
		vepDfFull = pd.DataFrame(columns = vepDfPart.columns)

	vepDfFull = pd.concat([vepDfFull, vepDfPart], ignore_index = True)


vepDfFull.to_csv("/net/data/aasubs/dbnsfp_catalog/pvs1/dbnsfp_stop_vep_output.csv")
