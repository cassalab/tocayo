import os


type_list = ["cv", "db"]


for df_type in type_list:

	for i in range(88):
		input_path = f"/net/data/aasubs/dbnsfp_catalog/vep/{i}_moderate_{df_type}_vep_input.vcf"
		output_path = f"/net/data/aasubs/dbnsfp_catalog/vep/{i}_moderate_{df_type}_vep_output.txt"

		job_path = f"/net/data/aasubs/dbnsfp_catalog/vep/jobs/j{i}_{df_type}_vep.sh"

		with open(job_path, "w") as f:

	 		command = f"""
			vep --cache --offline --dir /net/data/vep \
			--assembly GRCh38 \
			--tab --no_stats \
			-i {input_path} \
			-o {output_path} \
			--everything --force_overwrite --refseq --use_given_ref \
			--plugin CADD,/net/data/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz \
	 		"""

	 		f.write(command)

		f.close()

		os.system("qsub " + job_path)
