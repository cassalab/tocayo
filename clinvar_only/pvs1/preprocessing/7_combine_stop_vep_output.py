import pandas as pd


vepDfFull = pd.DataFrame()


for i in range(45):

	part_file_path = f"/net/data/aasubs/clinvar_only/pvs1/{i}_stop_vep_output.vcf"
	vepDfPart = pd.read_csv(part_file_path, sep = "\t", skiprows = 105)

	if i == 0:

		vepDfFull = pd.DataFrame(columns = vepDfPart.columns)

	vepDfFull = pd.concat([vepDfFull, vepDfPart], ignore_index = True)


vepDfFull.to_csv("/net/data/aasubs/clinvar_only/pvs1/stop_vep_output.csv")
