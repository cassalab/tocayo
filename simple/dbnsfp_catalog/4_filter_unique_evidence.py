import pandas as pd

evidenceList = ["strong", "moderate"]

for evidence in evidenceList:

	plpDfOriginal = pd.read_csv(f"dbnsfp_{evidence}_plp_analysis.csv")
	blbDfOriginal = pd.read_csv(f"dbnsfp_{evidence}_blb_analysis.csv")

	plpDf = plpDfOriginal[plpDfOriginal["cv_aa_half_name"].isin(blbDfOriginal["cv_aa_half_name"]) == False]
	blbDf = blbDfOriginal[blbDfOriginal["cv_aa_half_name"].isin(plpDfOriginal["cv_aa_half_name"]) == False]

	plpDf.drop(plpDf.filter(regex="Unname"),axis=1, inplace=True)
	blbDf.drop(blbDf.filter(regex="Unname"),axis=1, inplace=True)
	plpDf = plpDf.sort_values(by = ["cv_aa_sub_name"]).reset_index(drop = True)
	blbDf = blbDf.sort_values(by = ["cv_aa_sub_name"]).reset_index(drop = True)

	plpDf.to_csv(f"dbnsfp_{evidence}_plp_analysis.csv")
	blbDf.to_csv(f"dbnsfp_{evidence}_blb_analysis.csv")
