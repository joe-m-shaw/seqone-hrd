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

myriad_location_format_1 <- paste0(hrd_data_path, "myriad_reports_format_1/")

myriad_location_format_2 <- paste0(hrd_data_path, "myriad_reports_format_2/")

myriad_location_format_3 <- paste0(hrd_data_path, "myriad_reports_format_3/")

myriad_reports_format_1 <- list.files(myriad_location_format_1)

myriad_reports_format_2 <- list.files(myriad_location_format_2)

myriad_reports_format_3 <- list.files(myriad_location_format_3)

##################################################
# Functions
##################################################

read_myriad_pdf <- function(filepath, filename) {
  
  output <- pdftools::pdf_data(pdf = paste0(filepath, filename))
  
  return(output)
  
}

read_myriad_report_format_1 <- function(file) {
  
  myriad_report_data <- read_myriad_pdf(myriad_location_format_1, file)
  
  # Page 2: original Myriad report
  
  myriad_page2_text <- myriad_report_data[[2]][[6]]
  
  pg2_pathology_block <- paste0(myriad_page2_text[[82]], " ",
                                myriad_page2_text[[83]])
  
  pg2_igene_r_number <- myriad_page2_text[[76]]
  
  pg2_gi_score <- myriad_page2_text[[141]]
  
  pg2_brca_status <- myriad_page2_text[[159]]
  
  # Collate information together
  
  output <- data.frame("igene_r_number" = pg2_igene_r_number,
                       "gi_score" = pg2_gi_score,
                       "pathology_block" = pg2_pathology_block,
                       "brca_status" = pg2_brca_status,
                       "filename" = file)
  
  return(output)
  
}

read_myriad_report_format_2 <- function(file) {
  
  myriad_report_data <- read_myriad_pdf(myriad_location_format_2, file)
  
  myriad_page2_text <- myriad_report_data[[2]][[6]]
  
  myriad_page3_text <- myriad_report_data[[3]][[6]]
  
  pg2_igene_r_number <- myriad_page2_text[[75]]
  
  pg2_pathology_block <- paste0(myriad_page2_text[[81]], " ",
                                myriad_page2_text[[82]])
  
  pg3_gi_score <- myriad_page3_text[[77]]
  
  pg2_brca_status <- myriad_page2_text[[120]]
  
  output <- data.frame("igene_r_number" = pg2_igene_r_number,
                       "gi_score" = pg3_gi_score,
                       "pathology_block" = pg2_pathology_block,
                       "brca_status" = pg2_brca_status,
                       "filename" = file)
  
  return(output)
  
}

read_myriad_report_format_3 <- function(file) {
  
  myriad_report_data <- read_myriad_pdf(myriad_location_format_3, file)
  
  # Page 2: original Myriad report
  
  myriad_page2_text <- myriad_report_data[[2]][[6]]
  
  pg2_pathology_block <- paste0(myriad_page2_text[[83]], " ",
                                myriad_page2_text[[84]])
  
  pg2_igene_r_number <- myriad_page2_text[[77]]
  
  pg2_gi_score <- myriad_page2_text[[142]]
  
  pg2_brca_status <- myriad_page2_text[[161]]
  
  # Collate information together
  
  output <- data.frame("igene_r_number" = pg2_igene_r_number,
                       "gi_score" = pg2_gi_score,
                       "pathology_block" = pg2_pathology_block,
                       "brca_status" = pg2_brca_status,
                       "filename" = file)
  
  return(output)
  
}

##################################################
# Collate Report Information
##################################################

collated_info_format_1 <- data.frame()

for (i in myriad_reports_format_1) {
  
  tmp_output <- read_myriad_report_format_1(i)
  
  collated_info_format_1 <- rbind(collated_info_format_1, tmp_output)
  
  rm(tmp_output)
}


collated_info_format_2 <- data.frame()

for (i in myriad_reports_format_2) {
  
  tmp_output <- read_myriad_report_format_2(i)
  
  collated_info_format_2 <- rbind(collated_info_format_2, tmp_output)
  
  rm(tmp_output)
}

collated_info_format_3 <- data.frame()

for (i in myriad_reports_format_3) {
  
  tmp_output <- read_myriad_report_format_3(i)
  
  collated_info_format_3 <- rbind(collated_info_format_3, tmp_output)
  
  rm(tmp_output)
}

collated_myriad_results <- rbind(collated_info_format_1, collated_info_format_2,
                                 collated_info_format_3)

##################################################