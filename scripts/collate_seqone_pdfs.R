################################################################################
# Collating SeqOne HRD Report PDFs
# joseph.shaw3@nhs.net
################################################################################

##################################################
# Packages and Filepaths
##################################################

library(pdftools)
library(tidyverse)

source("scripts/hrd_filepaths.R")

##################################################
# Functions
##################################################

grep_seqone_text <- function(input_regex, page) {
  
  output <- sub(pattern = input_regex,
                 x = page,
                 replacement = "\\1")
  
  return(output)
  
}

read_seqone_report <- function(filepath, file) {
  
  seqone_report_text <- pdftools::pdf_text(pdf = paste0(filepath, file))
  
  page_1 <- seqone_report_text[[1]]
  
  hrd_score <- grep_seqone_text(".+CLASS\n.{78,83}((0.\\d{2})|(\\d{1})).+", page_1)
  
  hrd_status <- grep_seqone_text(".+SeqOne HRD Status1 : (\\D{8}).+", page_1)
  
  lga <- grep_seqone_text(".+LGA Status.{38}(\\d{1,2}).+", page_1)
  
  lpc <- grep_seqone_text(".+LPC Status.{38}(\\d{1,2}).+", page_1)
  
  ccne1 <- grep_seqone_text(".+CCNE1 Amplification.{29}((\\d{1}.\\d{2})|\\d{1}).+", page_1)
  
  rad51b <- grep_seqone_text(".+RAD51B Amplification.{28}(\\d{1}.\\d{1,2}).+", page_1)
  
  ncc <- grep_seqone_text(".+% of tumoral cells.{23,24}(\\d{2})%.+", page_1)
  
  coverage <- grep_seqone_text(".+Coverage.{33}(.{3,4})X.+", page_1)
  
  percent_mapping <- grep_seqone_text(".+% correct mapping.{24,25}(.{2,4})%.+", page_1)
  
  sample_id <- grep_seqone_text(".+Shallow sample ID.{15,18}\\D{2}\\d{6}_(.{8,26}).+", page_1)
  
  date <- grep_seqone_text(".+Date.{32}(\\D{3,9}.\\d{1,2}..\\d{4}).+", page_1)
  
  output <- data.frame(
    
    "hrd_score" = hrd_score,
    "hrd_status" = hrd_status,
    "lga" = lga,
    "lpc" = lpc,
    "ccne1" = ccne1,
    "rad51b" =rad51b,
    "ncc" = ncc,
    "coverage" = coverage,
    "percent_mapping" = percent_mapping,
    "sample_id" = sample_id,
    "date" = date,
    "filename" = i)
  
  return(output)
  
}

##################################################
# Locations and Files
##################################################

seqone_report_location <- paste0(hrd_data_path, "seqone_reports/")

seqone_reports <- list.files(seqone_report_location)

##################################################
# Collate Report Information
##################################################

collated_seqone_info <- data.frame()

for (i in seqone_reports) {
  
  tmp_output <- read_seqone_report(seqone_report_location, i)
  
  collated_seqone_info <- rbind(collated_seqone_info, tmp_output)
  
  rm(tmp_output)
  
}

##################################################