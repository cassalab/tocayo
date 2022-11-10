import pandas as pd


vepDfFull = pd.DataFrame()


for i in range(47):

	part_file_path = f"/net/data/aasubs/clinvar_only/vep/{i}_moderate_vep_output.txt"
	vepDfPart = pd.read_csv(part_file_path, sep = "\t", skiprows = 105)

	if i == 0:
		vepDfFull = pd.DataFrame(columns = vepDfPart.columns)

	vepDfFull = pd.concat([vepDfFull, vepDfPart], ignore_index = True)


vepDfFull.to_csv("/net/data/aasubs/clinvar_only/dup_files/moderate_vep_output.csv")
