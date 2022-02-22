import pandas as pd
import natsort
import sys
import os


evidence_list = ["strong", "moderate"]
type_list = ["cv", "db"]


for evidence in evidence_list:

	for df_type in type_list:

		if df_type == "cv":
			df = pd.read_csv(f"/net/data/aasubs/dbnsfp_catalog/{evidence}_{df_type}_potential_entries.csv", usecols = ["Chromosome", "Start", "simple_name", "ReferenceAlleleVCF", "AlternateAlleleVCF"], dtype = str)
		else:
			df = pd.read_csv(f"/net/data/aasubs/dbnsfp_catalog/{evidence}_{df_type}_potential_entries.csv", usecols = ["chr", "pos", "simple_name", "ref", "alt"], dtype = str)

		df["qual"] = "."
		df["filter"] = "."
		df["info"] = "."

		if df_type == "cv":
			df = df[["Chromosome", "Start", "simple_name", "ReferenceAlleleVCF", "AlternateAlleleVCF", "qual", "filter", "info"]]
		else:
			df = df[["chr", "pos", "simple_name", "ref", "alt", "qual", "filter", "info"]]

		index = natsort.index_natsorted(list(df["simple_name"]))

		for column in list(df.columns):
			df[column] = natsort.order_by_index(list(df[column]), index)

		headerList = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]

		headerText = [
		'##fileformat=VCFv4.2',
		'##FILTER=<ID=PASS,Description="All filters passed">',
		'##fileDate=20210817',
		'##reference=GRCh38/hg38',
		'##contig=<ID=1,length=248956422>',
		'##contig=<ID=2,length=242193529>',
		'##contig=<ID=3,length=198295559>',
		'##contig=<ID=4,length=190214555>',
		'##contig=<ID=5,length=181538259>',
		'##contig=<ID=6,length=170805979>',
		'##contig=<ID=7,length=159345973>',
		'##contig=<ID=8,length=145138636>',
		'##contig=<ID=9,length=138394717>',
		'##contig=<ID=10,length=133797422>',
		'##contig=<ID=11,length=135086622>',
		'##contig=<ID=12,length=133275309>',
		'##contig=<ID=13,length=114364328>',
		'##contig=<ID=14,length=107043718>',
		'##contig=<ID=15,length=101991189>',
		'##contig=<ID=16,length=90338345>',
		'##contig=<ID=17,length=83257441>',
		'##contig=<ID=18,length=80373285>',
		'##contig=<ID=19,length=58617616>',
		'##contig=<ID=20,length=64444167>',
		'##contig=<ID=21,length=46709983>',
		'##contig=<ID=22,length=50818468>',
		'##contig=<ID=X,length=156040895>',
		'##contig=<ID=Y,length=57227415>'
		]

		if evidence == "moderate":
			start = 0
			end = 70000

			for i in range(36):
				with open(f"/net/data/aasubs/dbnsfp_catalog/vep/{i}_{evidence}_{df_type}_vep_input.txt", "a") as f:
					
					for text in headerText:
						f.write(text)
						f.write("\n")

					f.write("\t".join(headerList))
					f.write("\n")

					if i == 35:
						for j in range(start, len(df)):
							f.write("\t".join(list(df.iloc[j])))
							f.write("\n")
					else:
						for j in range(start, end):
							f.write("\t".join(list(df.iloc[j])))
							f.write("\n")

				f.close()

				start += 70000
				end += 70000

				os.system(f"cp /net/data/aasubs/dbnsfp_catalog/vep/{i}_{evidence}_{df_type}_vep_input.txt /net/data/aasubs/dbnsfp_catalog/vep/{i}_{evidence}_{df_type}_vep_input.vcf")
				os.system(f"rm /net/data/aasubs/dbnsfp_catalog/vep/{i}_{evidence}_{df_type}_vep_input.txt")
		else:
			with open(f"/net/data/aasubs/dbnsfp_catalog/vep/{evidence}_{df_type}_vep_input.txt", "a") as f:
					
				for text in headerText:
					f.write(text)
					f.write("\n")

				f.write("\t".join(headerList))
				f.write("\n")

				for j in range(0, len(df)):
					f.write("\t".join(list(df.iloc[j])))
					f.write("\n")

				f.close()

				os.system(f"cp /net/data/aasubs/dbnsfp_catalog/vep/{evidence}_{df_type}_vep_input.txt /net/data/aasubs/dbnsfp_catalog/vep/{evidence}_{df_type}_vep_input.vcf")
				os.system(f"rm /net/data/aasubs/dbnsfp_catalog/vep/{evidence}_{df_type}_vep_input.txt")

