import pandas as pd
from natsort import natsorted

strong_cv_df = pd.read_csv("/net/data/aasubs/dbnsfp_catalog/strong_cv_potential_entries.csv", index_col = 0)
strong_db_df = pd.read_csv("/net/data/aasubs/dbnsfp_catalog/strong_db_potential_entries.csv", index_col = 0)

strong_cv_df.drop(strong_cv_df.filter(regex="Unname"), axis=1, inplace=True)
strong_db_df.drop(strong_db_df.filter(regex="Unname"), axis=1, inplace=True)

strong_cv_df["simple_name"] = strong_cv_df["Chromosome"].map(str) + " " + strong_cv_df["Start"].map(str) + " . " + strong_cv_df["ReferenceAlleleVCF"].map(str) + " " + strong_cv_df["AlternateAlleleVCF"].map(str) + " . . ."
strong_db_df["simple_name"] = strong_db_df["chr"].map(str) + " " + strong_db_df["pos"].map(str) + " . " +  strong_db_df["ref"].map(str) + " " + strong_db_df["alt"].map(str) + " . . ."

strong_cv_df["simple_name"] = natsorted(list(strong_cv_df["simple_name"]))
strong_db_df["simple_name"] = natsorted(list(strong_db_df["simple_name"]))

strong_cv_df["simple_name"].to_csv("/net/data/aasubs/dbnsfp_catalog/vep/strong_cv_vep_input.txt", sep = "\t", header = None, index = False)
strong_db_df["simple_name"].to_csv("/net/data/aasubs/dbnsfp_catalog/vep/strong_db_vep_input.txt", sep = "\t", header = None, index = False)
