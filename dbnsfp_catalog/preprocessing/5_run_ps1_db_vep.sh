vep --cache --offline --dir /net/data/vep --assembly GRCh38 --tab --no_stats -i /net/data/aasubs/update_11-06-22/dbnsfp_catalog/vep/strong_db_vep_input.vcf -o /net/data/aasubs/update_11-06-22/dbnsfp_catalog/vep/strong_db_vep_output.txt --everything --force_overwrite --refseq --use_given_ref --plugin CADD,/net/data/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz
