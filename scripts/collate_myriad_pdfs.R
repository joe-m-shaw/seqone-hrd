################################################################################
# Collating Myriad HRD Report PDFs
# joseph.shaw3@nhs.net
################################################################################

# Myriad HRD report PDFs have multiple formats
# Format 1: standard format
# Format 2: GI score is on page 2 of the Myriad report (page 3 of the whole
#           document)
# Format 3: "legal name" is used instead of first name and surname

##################################################
# Packages and Filepaths
##################################################

library(pdftools)
library(tidyverse)

source("scripts/hrd_filepaths.R")

##################################################
# Locations and Files
##################################################

myriad_reports_location <- paste0(hrd_data_path, "myriad_reports/")

myriad_report_files <- list.files(myriad_reports_location)

##################################################
# Functions
##################################################

read_myriad_pdf <- function(filepath, filename) {
  
  # Use pdf_text to read PDF as a single string per page
  
  output <- pdftools::pdf_text(pdf = paste0(filepath, filename))
  
  return(output)
  
}

read_myriad_report <- function(filepath, file) {
  
  myriad_report_text <- read_myriad_pdf(filepath, file)
  
  # iGene R Number
  
  r_number_regex <- ".+Patient ID:.{6}(R\\d{2}-\\w{4}).+"
  
  r_number <- sub(x = myriad_report_text[[2]],
                  pattern = r_number_regex,
                  replacement = "\\1")
  
  # NHS Number
  # Note: size of whitespace between "No:" and number can vary
  
  nhs_number_regex <- ".+NHS No:.{10,20}(\\d{3}.{1}\\d{3}.{1}\\d{4}).+"
  
  nhs_number <- sub(x = myriad_report_text[[1]],
                    pattern = nhs_number_regex,
                    replacement = "\\1")
  
  # Pathology Block
  # Note: some reports say "Block(s)" whilst others say "Specimen(s)"
  
  block_regex <- ".+(Block\\(s\\)|Specimen\\(s\\)) Analyzed:(.{5,20})\n\n.+"
  
  pathology_block <- sub(x = myriad_report_text[[2]],
                         pattern = block_regex,
                         replacement = "\\2")
  
  # Genomic Instability Score
  # Note: score can be 1 digit (i.e. "7") or 2 ("73")
  
  gi_score_regex <- ".+Patient Genomic Instability Score: (\\d{1,2}).+"
  
  gi_score_pg2 <- sub(x = myriad_report_text[[2]],
                  pattern = gi_score_regex,
                  replacement = "\\1")
  
  # Some reports have GI score on the third page
  
  gi_score_pg3 <- sub(x = myriad_report_text[[3]],
                      pattern = gi_score_regex,
                      replacement = "\\1")
  
  # Pick correct GI score
  
  gi_score <- ifelse(str_length(gi_score_pg2) %in% c(1,2),
                     gi_score_pg2,
                     ifelse(str_length(gi_score_pg3) %in% c(1,2),
                            gi_score_pg3, "NULL"))
  
  output <- data.frame("r_number" = r_number,
                       "nhs_number" = nhs_number,
                       "pathology_block" = pathology_block,
                       "gi_score" = gi_score)
  
  return(output)
  
}

##################################################
# Collate Report Information
##################################################

collated_myriad_info <- data.frame()

for (i in myriad_report_files) {
  
  tmp_output <- read_myriad_report(myriad_reports_location, i)
  
  collated_myriad_info <- rbind(collated_myriad_info, tmp_output)
  
  rm(tmp_output)
}

##################################################