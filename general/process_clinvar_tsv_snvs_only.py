import pandas as pd 


acceptable_benign_classifications = [
	"Benign",
	"Likely benign",
	"Benign/Likely benign"
]

acceptable_pathogenic_classifications = [
	"Pathogenic",
	"Likely pathogenic",
	"Pathogenic/Likely pathogenic"
]

acceptable_uncertain_classifications = [
	"Uncertain significance",
	"Conflicting interpretations of pathogenicity",
	"not provided"
]


def annototate_simply(annot):
	if annot in acceptable_benign_classifications:
		return "B/LB"
	elif annot in acceptable_pathogenic_classifications:
		return "P/LP"
	elif annot in acceptable_uncertain_classifications:
		return "VUS"
	else:
		return "Unrecognized"


def get_new_fields(df):
	aa_com = []
	aa_half = []
	refseqs = []
	stars = []
	for _, row in df.iterrows():
		status = row["ReviewStatus"]
		if status == "practice guideline":
			stars.append(4)
		elif status == "reviewed by expert panel":
			stars.append(3)
		elif status == "criteria provided, multiple submitters, no conflicts":
			stars.append(2)
		elif status == "criteria provided, conflicting interpretations" or status == "criteria provided, single submitter":
			stars.append(1)
		elif status == "no assertion for the individual variant" or status == "no assertion criteria provided" or status == "no assertion provided":
			stars.append(0)
		name = row["Name"]
		if "(p." in name:
			piece = name.split("(p.")[1]
			com = name.split(":")[1].split("(")[1][2:][:-1]
			half = com[:-3]
			refseq = name.split("(")[0]
			if "=" in piece or "?" in piece:
				aa_com.append(None)
				aa_half.append(None)
				refseqs.append(None)
			else:			
				aa_com.append(com)
				aa_half.append(half)
				refseqs.append(refseq)
		else:
			aa_com.append(None)
			aa_half.append(None)
			refseqs.append(None)

	df["aa_com"] = aa_com
	df["aa_half"] = aa_half
	df["refseq"] = refseqs
	df["star"] = stars
	return df


def parse_clinvar_df(fn = "/net/data/aasubs/variant_summary_11-06-22.txt"):
	df = pd.read_csv(fn, delimiter = "\t")
	df = df.loc[
		(~df["Start"].isna()) &
		(~df["ReferenceAlleleVCF"].isna()) &
		(~df["AlternateAlleleVCF"].isna())
	]

	df["simple_annot"] = df["ClinicalSignificance"].apply(annototate_simply)

	snv_df = df.loc[
		(df["Start"] == df["Stop"]) &
		(df["ReferenceAlleleVCF"] != df["AlternateAlleleVCF"]) &
		(df["simple_annot"] != "Unrecognized") &
		(df["Type"] == "single nucleotide variant") &
		(df["Assembly"] == "GRCh38")
	]

	snv_df["simple_name"] = (
		snv_df["Chromosome"].map(str) + "-" + 
		snv_df["Start"].map(str) + "-" + 
		snv_df["ReferenceAlleleVCF"].map(str) + "-" + 
		snv_df["AlternateAlleleVCF"].map(str)
	)

	return snv_df


if __name__ == "__main__":
	df = parse_clinvar_df().drop_duplicates(subset = ["simple_name"])
	df = get_new_fields(df).dropna(subset = ["aa_com"])
	df["gene_name"] = df["Name"].str.split("(", 1, expand = True)[1].str.split(")", 1, expand = True)[0]
	df["aa_sub_name"] = df["gene_name"] + "-" + df["aa_com"]
	df["aa_half_name"] = df["gene_name"] + "-" + df["aa_half"]
	df["Chromosome"] = df["Chromosome"].astype(str)
	df.drop(df.filter(regex="Unname"), axis=1, inplace=True)
	df.sort_values(by = ["aa_sub_name"]).reset_index(drop = True).to_csv("/net/data/aasubs/parsed_clinvar_11-06-22.csv")
