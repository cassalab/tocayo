import pandas as pd
import sys


type_list = ["cv", "db"]


for df_type in type_list:

	vepDfFull = pd.DataFrame()

	for i in range(88):
		part_file_path = f"/net/data/aasubs/dbnsfp_catalog/vep/{i}_moderate_{df_type}_vep_output.txt"

		vepDfPart = pd.read_csv(part_file_path, sep = "\t", skiprows = 105)

		if i == 0:
			vepDfFull = pd.DataFrame(columns = vepDfPart.columns)

		vepDfFull = pd.concat([vepDfFull, vepDfPart], ignore_index = True)


	vepDfFull.to_csv(f"/net/data/aasubs/dbnsfp_catalog/vep/moderate_{df_type}_vep_output.csv")
