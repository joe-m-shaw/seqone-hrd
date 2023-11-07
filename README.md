# Homologous Recombination Deficiency Validation

This repo contains the analysis script for the validation of Homologous Recombination Deficiency (HRD) testing using SeqOne shallow whole genome sequencing at the North West Genomic Laboratory Hub (Manchester).

## Data

Data were downloaded from the SeqOne website (https://seqone.com) as PDFs, and the pdftools package was used to convert them into dataframes. SeqOne do not currently offer the option of downloading information in a tabulated form.

## Scripts

The hrd_analysis.R script includes all analysis for the validation and should run with the associated data (the associated data is not included in this repo). The separate dlms_queries.R script was used to extract sample information from the DNA Laboratory Management System (DLMS) via an ODBC connection, which is set up for PC38698. The tables were then saved as csvs. The dlms_queries.R script will not run on a computer not set up with an ODBC DLMS connection, and is not required to rerun the hrd_analysis.R script.

## SomaHRD Versions

During the validation, versions 1.1. and 1.2 of the SomaHRD pipeline were examined. Full details of the differences between the versions are provided in the validation document. There are some differences in report formatting between versions, and the "robustness" and "low tumour fraction" metrics are only provided for version 1.2.
