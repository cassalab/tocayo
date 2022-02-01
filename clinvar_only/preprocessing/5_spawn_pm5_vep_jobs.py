import os

files_left = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"]

for x in files_left:
	input_path = f"/net/data/aasubs/clinvar_only/vep/{x}_moderate_vep_input.txt"
	output_path = f"/net/data/aasubs/clinvar_only/vep/{x}_moderate_vep_output.txt"
	plugin_command = "--plugin CADD,/net/data/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz"
	assembly = "GRCh38"

	job_path = "/net/data/aasubs/clinvar_only/vep/j" + x + "_vep.sh"
	with open(job_path, "w") as f:
 		command = f"""
 		vep --cache --dir /net/data/vep --assembly {assembly} --tab --no_stats -i {input_path} -o {output_path} --everything --refseq {plugin_command}
 		"""
 		f.write(command)
	f.close()

	os.system("qsub " + job_path)
