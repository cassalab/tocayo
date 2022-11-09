import os

for i in range(47):
	input_path = f"/net/data/aasubs/clinvar_only/vep/{i}_moderate_vep_input.txt"
	output_path = f"/net/data/aasubs/clinvar_only/vep/{i}_moderate_vep_output.txt"
	plugin_command = "--plugin CADD,/net/data/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz"
	assembly = "GRCh38"

	job_path = "/net/data/aasubs/clinvar_only/vep/j" + str(i) + "_vep.sh"
	with open(job_path, "w") as f:
 		command = f"""
 		vep --cache --offline --dir /net/data/vep --assembly {assembly} --tab --no_stats -i {input_path} -o {output_path} --everything --refseq --use_given_ref {plugin_command}
 		"""
 		f.write(command)
	f.close()

	os.system("qsub " + job_path)
