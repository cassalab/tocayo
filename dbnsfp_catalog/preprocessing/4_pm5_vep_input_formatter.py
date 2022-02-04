import pandas as pd
from natsort import natsorted

moderate_cv_df = pd.read_csv("/net/data/aasubs/dbnsfp_catalog/moderate_cv_potential_entries.csv", index_col = 0)
moderate_db_df = pd.read_csv("/net/data/aasubs/dbnsfp_catalog/moderate_db_potential_entries.csv", index_col = 0)

moderate_cv_df.drop(moderate_cv_df.filter(regex="Unname"), axis=1, inplace=True)
moderate_db_df.drop(moderate_db_df.filter(regex="Unname"), axis=1, inplace=True)

moderate_cv_df["simple_name"] = moderate_cv_df["Chromosome"].map(str) + " " + moderate_cv_df["Start"].map(str) + " . " + moderate_cv_df["ReferenceAlleleVCF"].map(str) + " " + moderate_cv_df["AlternateAlleleVCF"].map(str) + " . . ."
moderate_db_df["simple_name"] = moderate_db_df["chr"].map(str) + " " + moderate_db_df["pos"].map(str) + " . " +  moderate_db_df["ref"].map(str) + " " + moderate_db_df["alt"].map(str) + " . . ."

moderate_cv_df["simple_name"] = natsorted(list(moderate_cv_df["simple_name"]))
moderate_db_df["simple_name"] = natsorted(list(moderate_db_df["simple_name"]))

moderate_cv_df["simple_name"].to_csv("/net/data/aasubs/dbnsfp_catalog/vep/moderate_cv_vep_input.txt", sep = "\t", header = None, index = False)
moderate_db_df["simple_name"].to_csv("/net/data/aasubs/dbnsfp_catalog/vep/moderate_db_vep_input.txt", sep = "\t", header = None, index = False)
