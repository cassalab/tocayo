# Tocayo

Tocayo is a pipeline for analyzing sequence variants which have evidence of pathogenicity according to the PS1, PM5, and PVS1 guidelines specified by the ACMG/AMP. Specifically, the pipeline is meant to analyze VUS in ClinVar and other unclassified nsSNVs.

---
## Downloads (recommended)

Pre-computed lists of variants (VUS within ClinVar and unclassified variants exome-wide) with PS1, PM5, and PVS1 evidence are accessible through the following links [created using ClinVar v12/19/2021 and dbNSFP v4.2a]:
* Excel Version - https://figshare.com/s/a8ed691e7985afdbb0a1
* CSV Version - https://figshare.com/s/b29d7d45d62bff0b3bc9

All lists (Excel and CSV) are updated every 6 months with the latest version of the ClinVar variant summary TSV file. Usage of these lists is recommended for all purposes unless having the most current version of ClinVar is of great importance. 

---
## Pipeline Requirements
To run the pipeline on an updated version of ClinVar, the following files and tools are required:
* ClinVar variant_summary.txt file, available at ClinVar's Downloads/FTP site (https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/)
* dbNSFP .sqlite file with relavent annotations for all ~84M nsSNVs, available at https://figshare.com/s/daadeb328a109f1f5f65
* Ensembl VEP v108 downloaded locally with the CADD and LOFTEE plugins

_Note:_
_Using dbNSFP files with updated annotations is not necessary as only fundamental variant data is used within the pipeline._


Unzip the dbNSFP .sqlite.gz input file before use.

```
gunzip dbnsfp_subset.sqlite.gz
```

---
## Usage

The files above can be run sequentially to obtain lists of variants with PS1, PM5, and PVS1 evidence. Contents in /clinvar_only/ and /dbnsfp_catalog/ are for obtaining VUS variants within ClinVar and unclassified variants exome-wide with ACMG/AMP evidence, respectively. ClinVar variant summary .txt files should be run through the processor in /general/, then used to run subsequent files in /preprocessing/ and /analysis/ folders in this specific order. Replace all file paths with appropriate local file paths; outputs from one file are often used in subsequent files. Using a similar file path structure and similar file names is suggested.

_Note:_ Ensemble VEP is needed to retrieve annotations (with refseq option enabled and .txt output). The CADD and LOFTEE plugins for VEP are also used at various points. More information on these plugins can be found here: https://ensembl.org/info/docs/tools/vep/script/vep_plugins.html
