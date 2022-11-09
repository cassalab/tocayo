import pandas as pd
import os

df = pd.read_csv("/net/data/aasubs/update_11-06-22/clinvar_only/dup_files/moderate_vep_input.txt", sep = "\t", index_col = False, header = None)

start = 0
end = 5000

for i in range(47):
	with open(f"/net/data/aasubs/update_11-06-22/clinvar_only/vep/{i}_moderate_vep_input.txt", "a") as f:
		if i == 46:
			for j in range(start, len(df)):
				f.write(str(df.iloc[j][0]))
				f.write("\n")
		else:
			for j in range(start, end):
				f.write(str(df.iloc[j][0]))
				f.write("\n")
	f.close()

	start += 5000
	end += 5000
