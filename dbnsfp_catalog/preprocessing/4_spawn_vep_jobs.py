import os


files_left = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35"]

type_list = ["cv", "db"]


for df_type in type_list:

	for x in files_left:
		input_path = f"/net/data/aasubs/dbnsfp_catalog/vep/{x}_moderate_{df_type}_vep_input.vcf"
		output_path = f"/net/data/aasubs/dbnsfp_catalog/vep/{x}_moderate_{df_type}_vep_output.vcf"

		job_path = f"/net/data/aasubs/dbnsfp_catalog/vep/jobs/j{x}_{df_type}_vep.sh"

		with open(job_path, "w") as f:

	 		command = f"""
			vep --cache --offline --dir /net/data/vep \
			--assembly GRCh38 \
			--fasta /net/data/vep/homo_sapiens/105_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa \
			--tab --no_stats \
			-i {input_path} \
			-o {output_path} \
			--everything --force_overwrite --refseq \
			--plugin CADD,/net/data/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz \
	 		"""

	 		f.write(command)

		f.close()

		os.system("qsub " + job_path)
