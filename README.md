# Homologous Recombination Deficiency Testing 

This repo contains scripts for validating and monitoring Homologous Recombination Deficiency (HRD) testing using SeqOne shallow whole genome sequencing at the North West Genomic Laboratory Hub (GLH) in Manchester.

## Project structure

### data

HRD results from the SeqOne SomaHRD pipeline are downloaded as csvs or PDFs the SeqOne website (https://seqone.com). All results are saved internally at the North West GLH. **No data should be available in this Github repository.**

The filepath for the data is specified in the config.yml file and called in R using `config::get()`.

### functions

Functions for reading and manipulating HRD results are saved in the functions folder.

### scripts

The scripts folder contains scripts that process, collate and reformat HRD results. 

### vignettes

Final reports for validations (DOC prefix) or incident investigations (INC prefix) are prepared as Quarto markdown documents and saved in the vignettes folder. The documents are named with the document identifiers from the Manchester GLH quality management system, QPulse.
