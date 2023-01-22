import pandas as pd


vepDfFull = pd.DataFrame()


for i in range(102):
	part_file_path = f"/net/data/aasubs/dbnsfp_catalog/pvs1/{i}_dbn_splice_vep_output.vcf"
	vepDfPart = pd.read_csv(part_file_path, sep = "\t", skiprows = 109)

	if i == 0:
		vepDfFull = pd.DataFrame(columns = vepDfPart.columns)

	vepDfFull = pd.concat([vepDfFull, vepDfPart], ignore_index = True)


vepDfFull.to_csv("/net/data/aasubs/dbnsfp_catalog/pvs1/dbnsfp_splice_vep_output.csv")
