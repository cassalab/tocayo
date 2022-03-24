import pandas as pd

evidenceList = ["strong", "moderate"]
classificationList = ["P/LP", "B/LB"]
classificationList_v2 = ["plp", "blb"]

for evidence in evidenceList:

	for i in range(len(classificationList)):

		clas = classificationList[i]
		clas_v2 = classificationList_v2[i]
		
		cvDf = pd.read_csv(f"/net/data/aasubs/dbnsfp_catalog/{evidence}_cv_potential_entries.csv")
		dbDf = pd.read_csv(f"/net/data/aasubs/dbnsfp_catalog/{evidence}_db_potential_entries.csv")
		
		dbDf["aa_sub_name"] = cvDf["aa_sub_name"]
		
		cvDf.drop(cvDf.filter(regex="Unname"), axis=1, inplace=True)
		dbDf.drop(dbDf.filter(regex="Unname"), axis=1, inplace=True)

		cv_columns = []
		for column in list(cvDf.columns):
			cv_columns.append("cv_" + column)

		db_columns = []
		for column in list(dbDf.columns):
			db_columns.append("db_" + column)

		cvDf.columns = cv_columns
		dbDf.columns = db_columns

		df = pd.concat([cvDf, dbDf], axis = 1)
		df = df.drop_duplicates(subset = ["cv_simple_name", "db_simple_name"])

		df = df.sort_values(by = ["cv_aa_sub_name"]).reset_index(drop = True)

		df = df[df["cv_simple_annot"] == clas]

		df.to_csv(f"/net/data/aasubs/dbnsfp_catalog/classified/dbnsfp_{evidence}_{clas_v2}_analysis.csv")
