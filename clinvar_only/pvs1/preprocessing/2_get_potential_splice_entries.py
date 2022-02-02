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


def annotate_star(status):
	if status == "practice guideline":
		return 4
	elif status == "reviewed by expert panel":
		return 3
	elif status == "criteria provided, multiple submitters, no conflicts":
		return 2
	elif status == "criteria provided, conflicting interpretations" or status == "criteria provided, single submitter":
		return 1
	elif status == "no assertion for the individual variant" or status == "no assertion criteria provided" or status == "no assertion provided":
		return 0
	else:
		return -1


if __name__ == '__main__':

	df = pd.read_csv("/net/data/aasubs/variant_summary.txt", sep = "\t")

	df = df.loc[
		(~df["Start"].isna()) &
		(~df["ReferenceAlleleVCF"].isna()) &
		(~df["AlternateAlleleVCF"].isna())
	]

	df["simple_annot"] = df["ClinicalSignificance"].apply(annototate_simply)
	df["star"] = df["ReviewStatus"].apply(annotate_star)

	df = df.loc[
		(df["Start"] == df["Stop"]) &
		(df["ReferenceAlleleVCF"] != df["AlternateAlleleVCF"]) &
		(df["simple_annot"] != "Unrecognized") &
		(df["Type"] == "single nucleotide variant") &
		(df["Assembly"] == "GRCh38") &
		(df["Name"].str.contains("=") == False)
	]

	df["simple_name"] = (
		df["Chromosome"].map(str) + "-" + 
		df["Start"].map(str) + "-" + 
		df["ReferenceAlleleVCF"].map(str) + "-" + 
		df["AlternateAlleleVCF"].map(str)
	)

	df["refseq"] = df["Name"].str.split("(").get(0)
	df["gene_name"] = df["Name"].str.split("(").str.get(1).str.split(")").str.get(0)
	df["Chromosome"] = df["Chromosome"].astype(str)
	df.drop(df.filter(regex="Unname"), axis=1, inplace=True)

	df = df[(df["Name"].str[-5:-3] == "-1") | (df["Name"].str[-5:-3] == "-2") | (df["Name"].str[-5:-3] == "+1") | (df["Name"].str[-5:-3] == "+2")]
	
	df.to_csv("/net/data/aasubs/clinvar_only/pvs1/splice_entries.csv")

