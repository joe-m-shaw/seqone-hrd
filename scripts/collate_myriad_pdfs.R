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
# Functions
##################################################

read_myriad_report <- function(filepath, file) {
  
  # Use pdf_text to read PDF as a single string per page
  
  myriad_report_text <- pdftools::pdf_text(pdf = paste0(filepath, file))
  
  page1 <- myriad_report_text[[1]]
  
  page2 <- myriad_report_text[[2]]
  
  page3 <- myriad_report_text[[3]]
  
  # iGene R Number
  
  myriad_r_number <- sub(x = page2,
                  pattern = ".+Patient ID:.{6}(R\\d{2}-\\w{4}).+",
                  replacement = "\\1")
  
  # NHS Number
  # Note: size of whitespace between "No:" and number can vary
  
  nhs_number <- sub(x = page1,
                    pattern = ".+NHS No:.{10,20}(\\d{3}.{1}\\d{3}.{1}\\d{4}).+",
                    replacement = "\\1")
  
  nhs_number_mod <- as.numeric(gsub(pattern = "(\\D)", "",
                                    nhs_number))
  
  
  # Patient name
  myriad_patient_name <- sub(x = page1,
                              pattern = ".+Patient Name:\\s{10,15}(\\D{5,25})\n.+",
                              replacement = "\\1")
  
  # Patient date of birth
  myriad_dob <- sub(x = page1,
                    pattern = ".+Date of Birth:\\s{10,15}(\\d{2}/\\d{2}/\\d{4})\n.+",
                    replacement = "\\1")
  
  # Pathology Block
  # Note: some reports say "Block(s)" whilst others say "Specimen(s)"
  
  myriad_pathology_block_pg2 <- gsub("[\n]", "", sub(x = page2,
                         pattern = ".+(Block\\(s\\)|Specimen\\(s\\)) Analyzed:(.{5,20})\n\n.+",
                         replacement = "\\2"))
  
  myriad_pathology_block_pg1 <- gsub("[\n]", "", sub(x = page1,
                                    pattern = ".+Pathology No:\\s{10,30}(.{5,25})(\\n\\n\\n\\n\\n|\\n).+",
                                    replacement = "\\1"))
  
  # Genomic Instability Score
  # Note: score can be 1 digit (i.e. "7") or 2 ("73")
  
  gi_score_regex <- ".+Patient Genomic Instability Score: (\\d{1,2}).+"
  
  gi_score_pg2 <- as.numeric(sub(x = page2,
                  pattern = gi_score_regex,
                  replacement = "\\1"))
  
  # Some reports have GI score on the third page
  
  gi_score_pg3 <- as.numeric(sub(x = page3,
                      pattern = gi_score_regex,
                      replacement = "\\1"))
  
  # Pick correct GI score
  
  myriad_gi_score <- ifelse(str_length(gi_score_pg2) %in% c(1,2),
                     gi_score_pg2,
                     ifelse(str_length(gi_score_pg3) %in% c(1,2),
                            gi_score_pg3, "NULL"))
  
  myriad_hrd_status <- sub(x = page2,
                           pattern = ".+Myriad HRD Status:.(\\D{8}).+",
                           replacement = "\\1")
  
  myriad_brca_status <- sub(x = page2,
                          pattern = ".+Tumor Mutation BRCA1/BRCA2 Status:.(\\D{8}).+",
                          replacement = "\\1")
  
  
  output <- data.frame("myriad_r_number" = myriad_r_number,
                       "myriad_patient_name" = myriad_patient_name,
                       "myriad_dob" = myriad_dob,
                       "nhs_number" = nhs_number_mod,
                       "myriad_pathology_block_pg1" = myriad_pathology_block_pg1,
                       "myriad_pathology_block_pg2" = myriad_pathology_block_pg2,
                       "myriad_gi_score" = myriad_gi_score,
                       "myriad_hrd_status" = myriad_hrd_status,
                       "myriad_brca_status" = myriad_brca_status,
                       "myriad_filename" = file)
  
  return(output)
  
}

##################################################
# Locations and Files
##################################################

myriad_reports_location <- paste0(hrd_data_path, "myriad_reports/")

myriad_report_files <- list.files(myriad_reports_location)

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
# Checks
##################################################

# Check for NA values
stopifnot(collated_myriad_info[is.na(collated_myriad_info)] == 0)

# Check all files included
stopifnot(setdiff(myriad_report_files, collated_myriad_info$myriad_filename) == 0)

# Check GI score range is as expected
stopifnot(
  
  max(collated_myriad_info$myriad_gi_score) <= 100,
  
  min(collated_myriad_info$myriad_gi_score) >= 0
  
)

# Check HRD status is dichotomous
stopifnot(setdiff(unique(collated_myriad_info$myriad_hrd_status),
                  c("POSITIVE", "NEGATIVE")) == 0)

##################################################