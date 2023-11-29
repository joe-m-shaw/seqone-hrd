# Homologous Recombination Deficiency Validation

This repo contains the analysis script for the validation of Homologous Recombination Deficiency (HRD) testing using SeqOne shallow whole genome sequencing at the North West Genomic Laboratory Hub (Manchester).

## Data

Reports were downloaded from the SeqOne website (https://seqone.com) as PDFs, and the pdftools package was used to convert them into dataframes. SeqOne do not currently offer the option of downloading information in a tabulated form.

The separate dlms_queries.R script was used to extract sample information from the DNA Laboratory Management System (DLMS) via an ODBC connection, which is set up for PC38698. The tables were then saved as csvs. The dlms_queries.R script will not run on a computer not set up with an ODBC DLMS connection, and is not required to rerun the hrd_analysis.R script. 

All the data used in this validation is saved in the "data" folder at this interal location: S:\central shared\Genetics\Mol_Shared\Development.Team\SeqOne Homologous Recombination Deficiency Validation\HRD R script files\data

## Script

The hrd_analysis.R script includes all analysis for the validation and uses local filepaths.

## Running the Analysis

To recreate the analysis performed for the validation, copy the data folder into the project folder. You will also need to create "outputs" and "plots" folders.

The hrd_analysis.R analysis script should then produce the plots and tables included in the final validation document.

## SomaHRD Versions

During the validation, versions 1.1. and 1.2 of the SomaHRD pipeline were examined. Full details of the differences between the versions are provided in the validation document. There are some differences in report formatting between versions, and the "robustness" and "low tumour fraction" metrics are only provided for version 1.2.
