################################################################################
# SeqOne HRD Analysis
# joseph.shaw3@nhs.net
################################################################################

rm(list=ls())

##################################################
# Packages
##################################################

library(pdftools)
library(tidyverse)

##################################################
# Read PDFs
##################################################

location <- "~/homologous_recombination_deficiency/inputs/"

report <- pdftools::pdf_text(pdf = paste0(location, "seqone-hrd-example-report-positive-brca", ".pdf")) %>% 
  strsplit(split = "\n")

pdftools::pdf_data(pdf = paste0(location, "seqone-hrd-example-report-positive-brca", ".pdf"))

##################################################