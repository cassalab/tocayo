import os

files_left = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18"]

for x in files_left:
	input_path = f"/net/data/aasubs/clinvar_only/pvs1/{x}_splice_vep_input.vcf"
	output_path = f"/net/data/aasubs/clinvar_only/pvs1/{x}_splice_vep_output.vcf"

	job_path = "/net/data/aasubs/clinvar_only/pvs1/j" + x + "_splice_vep_lof.sh"

	with open(job_path, "w") as f:

 		command = f"""
 		#!/bin/bash
		eval $(perl -Mlocal::lib=--deactivate-all)
		eval "$(/net/apps/x86_64/miniconda3/bin/conda shell.bash hook)"
		vep --cache --offline --dir /net/data/vep \
		--assembly GRCh38 \
		--fasta /net/data/vep/homo_sapiens/105_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa \
		--tab --no_stats \
		-i {input_path} \
		-o {output_path} \
		--everything --force_overwrite --refseq \
		--plugin CADD,/net/data/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz \
		--plugin LoF,loftee_path:/net/data/vep/loftee-grch38,human_ancestor_fa:/net/data/vep/loftee-grch38/human_ancestor.fa.gz,conservation_file:/net/data/vep/loftee-grch38/loftee.sql,gerp_bigwig:/net/data/vep/loftee-grch38/gerp_conservation_scores.homo_sapiens.GRCh38.bw \
		--dir_plugins /net/data/vep/loftee-grch38
 		"""

 		f.write(command)

	f.close()

	os.system("qsub " + job_path)
