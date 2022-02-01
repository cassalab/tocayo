import pandas as pd
from natsort import natsorted

strong_df = pd.read_csv("/net/data/aasubs/clinvar_only/dup_files/strong_all_dup_subs.csv", index_col = 0)
moderate_df = pd.read_csv("/net/data/aasubs/clinvar_only/dup_files/moderate_all_dup_subs.csv", index_col = 0)

strong_df.drop(strong_df.filter(regex="Unname"), axis=1, inplace=True)
moderate_df.drop(moderate_df.filter(regex="Unname"), axis=1, inplace=True)

strong_df["simple_name"] = strong_df["Chromosome"].map(str) + " " + strong_df["Start"].map(str) + " . " + strong_df["ReferenceAlleleVCF"].map(str) + " " + strong_df["AlternateAlleleVCF"].map(str) + " . . ."
moderate_df["simple_name"] = moderate_df["Chromosome"].map(str) + " " + moderate_df["Start"].map(str) + " . " +  moderate_df["ReferenceAlleleVCF"].map(str) + " " + moderate_df["AlternateAlleleVCF"].map(str) + " . . ."

strong_df["simple_name"] = natsorted(list(strong_df["simple_name"]))
moderate_df["simple_name"] = natsorted(list(moderate_df["simple_name"]))

strong_df["simple_name"].to_csv("/net/data/aasubs/clinvar_only/dup_files/strong_vep_input.txt", sep = "\t", header = None, index = False)
moderate_df["simple_name"].to_csv("/net/data/aasubs/clinvar_only/dup_files/moderate_vep_input.txt", sep = "\t", header = None, index = False)
