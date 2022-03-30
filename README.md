# Tocayo

## Usage

The files above can be run sequentially to obtain lists of variants with PS1, PM5, and PVS1 evidence. Contents in /clinvar_only/ and /dbnsfp_catalog/ are for obtaining VUS variants within ClinVar and unclassified variants exome-wide with ACMG/AMP evidence, respectively. ClinVar variant summary .txt files should be run through the processor in /general/, then used to run subsequent files in /preprocessing/ and /analysis/ folders in this specific order. Replace all file paths with appropriate local file paths; outputs from one file are often used in subsequent files. Using a similar file path structure and similar file names are suggested.

A far simpler and smaller set of files to obtain similar lists of variants are available in /simple/. However, amino acid substitutions are not verified to be the same or in the same codon for PS1 and PM5 evidences, consequences are not checked for all evidence tiers, and unclassified variants from dbNSFP cannot be matched by transcript for PS1 and PM5 evidence. It is recommended to use this set of files when ease of running is a priority and when additional annotations are not needed, and it is recommended to use the main set of files when having the most accurate data with all annotations is a priority.

_Note:_ When using the main set of files, Ensemble VEP is needed to retrieve annotations (with refseq option enabled and .txt output). CADD and LOFTEE plugins are used at various points, but these aspects are inessential and can be removed from respective files. Nothing outside of input files is needed to run the simple set of files.

---
## Downloads

Pre-computed lists of variants (VUS within ClinVar and unclassified variants exome-wide) with PS1, PM5, and PVS1 evidence are accessible through the following links [created using ClinVar v12/19/2021 and dbNSFP v4.2a]:
* Excel Version - https://figshare.com/s/a8ed691e7985afdbb0a1
* CSV Version - https://figshare.com/s/b29d7d45d62bff0b3bc9

ClinVar and dbNSFP files which can be inputted into this pipeline are accessible through the following link:
* https://figshare.com/s/daadeb328a109f1f5f65

_Note:_ It is recommended to use the most recent variant_summary.txt file from ClinVar's Downloads/FTP site (https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/).
Using dbNSFP files with updated annotations is not necessary as only fundamental variant data from the databse is used within the pipeline.
