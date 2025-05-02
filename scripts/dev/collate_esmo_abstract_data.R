# ESMO abstract data collation

# May 2025: Dr Rob Morgan is preparing an abstract for the ESMO conference which 
# will include data from the Manchester SeqOne service.

library(tidyverse)

source(here::here("functions/hrd_functions.R"))

esmo_folderpath <- paste0(config::get("data_filepath"),
                          "live_service/ESMO_paper/")

# DNA Database connection -------------------------------------------------

message("Connecting to DNA database")

dbi_con <- DBI::dbConnect(
  drv = odbc::odbc(),
  dsn = "moldb")

sample_tbl <- tbl(dbi_con, dbplyr::in_catalog(catalog = "MolecularDB",
                                              schema = "dbo",
                                              table = "Samples")) |> 
  janitor::clean_names()

dna_db_worksheets <- tbl(dbi_con, dbplyr::in_catalog(catalog = "MolecularDB",
                                                     schema = "dbo",
                                                     table = "PCR_New"))|> 
  janitor::clean_names()

# Get list of SeqOne HRD worksheets ---------------------------------------

message("Retrieving SeqOne worksheet list")

all_worksheet_df <- dna_db_worksheets |> 
  select(pcrid, description) |> 
  collect()

seqone_worksheet_df <- all_worksheet_df |> 
  filter(grepl(pattern = "seqone|seq_one|seq\\sone",
               x = description,
               ignore.case = TRUE)) |> 
  mutate(worksheet = paste0("WS", pcrid))

seqone_worksheet_list <- list(seqone_worksheet_df$worksheet)

stopifnot(length(seqone_worksheet_list[[1]]) > 1)

# Find csvs in worksheet folders ------------------------------------------

message("Locating SeqOne csv files")

find_seqone_csvs <- function(worksheet) {
  
  data_path <- "S:/central shared/Genetics/Repository/WorksheetAnalysedData/"
  
  folder_path <- str_c(data_path, worksheet, "/")
  
  output <- list.files(folder_path,
                       full.names = TRUE,
                       recursive = TRUE,
                       pattern = "hrd-results.*csv")
  
}

seqone_worksheet_filepaths <- seqone_worksheet_list |> 
  map(\(seqone_worksheet_list) find_seqone_csvs(seqone_worksheet_list)) |> 
  flatten()

# Copy to working directory -----------------------------------------------

message("Copying csv files to directory")

file.copy(from = seqone_worksheet_filepaths,
          to = paste0(esmo_folderpath, "raw/"))

# Collate csv files -------------------------------------------------------

message("Collating files")

files_to_collate <- list.files(paste0(esmo_folderpath, "raw/"),
                     full.names = TRUE,
                     recursive = FALSE,
                     pattern = "hrd-results.*csv")

collated_seqone_csv_data <- files_to_collate |> 
  map(\(files_to_collate) read_seqone_csv(file = files_to_collate)) |> 
  list_rbind() |> 
  # Remove two samples with poor quality and absent scores
  filter(!is.na(score))

collated_seqone_csv_data_mod <- collated_seqone_csv_data |> 
  mutate(date = parse_date_time(x = analysis_date, 
                       orders = c("dmy", "ymd")),
         worksheet = str_extract(string = sample,
                                 pattern = "(WS[0-9]{6})_(\\d{8})",
                                 group = 1),
         labno = str_extract(string = sample,
                             pattern = "(WS[0-9]{6})_(\\d{8})",
                             group = 2),
         status = factor(status,
                         levels = c("Positive", "Negative",
                                    "Non-conclusive")))

# Check data --------------------------------------------------------------

message("Checking data")

stopifnot(anyNA.data.frame(collated_seqone_csv_data_mod |> 
                             select(-c(brca_status, brca_mutation))) == FALSE)

if((min(collated_seqone_csv_data_mod$LGA) < 0) |
   (max(collated_seqone_csv_data_mod$LGA) > 100)){
  stop("LGA should be in range 0-100")
} else {
  message("LGA check passed")
}

if((min(collated_seqone_csv_data_mod$LPC) < 0)|
   (max(collated_seqone_csv_data_mod$LPC) > 100)){
     stop("LPC should be in range 0-100")
} else {
  message("LPC check passed")
}

if((min(collated_seqone_csv_data_mod$score) < 0)|
   (max(collated_seqone_csv_data_mod$score) > 1)){
  stop("SeqOne score should be between -10")
} else{
  message("SeqOne score check passed")
}

stopifnot(levels(collated_seqone_csv_data_mod$status) == 
            c("Positive", "Negative", "Non-conclusive"))

if((min(collated_seqone_csv_data_mod$ccne1_cn) < 0) |
   max(collated_seqone_csv_data_mod$ccne1_cn) > 100){
  stop("Check CCNE1 copy number")
} else {
  message("CCNE1 copy number check passed")
}

if((min(collated_seqone_csv_data_mod$rad51b_cn) < 0) |
    (max(collated_seqone_csv_data_mod$rad51b_cn) > 100)){
  stop("Check RAD51B copy number")
} else {
  message("RAD51B copy number check passed")
}

if((min(collated_seqone_csv_data_mod$coverage) < 0) |
   (max(collated_seqone_csv_data_mod$coverage) > 7)){
  stop("Check coverage")
} else {
  message("Coverage check passed")
}

if((min(collated_seqone_csv_data_mod$pct_tum_cell) < 0) |
   (max(collated_seqone_csv_data_mod$pct_tum_cell) > 1)){
  stop("Check percentage tumour cells")
} else {
  message("Percentage tumour cells check passed")
}

if((min(collated_seqone_csv_data_mod$gi_confidence) < 0) |
   (max(collated_seqone_csv_data_mod$gi_confidence) > 1)){
  stop("Check GI confidence")
} else {
  message("GI confidence check passed")
}

message("Data check complete")

# Add patient identifiers -------------------------------------------------

message("Adding patient identifiers")

seqone_labnos <- unique(collated_seqone_csv_data_mod$labno)

seqone_sample_info <- sample_tbl |> 
  filter(labno %in% seqone_labnos) |> 
  select(labno, firstname, surname, nhsno, dob, comments) |> 
  collect()

seqone_csv_data_with_ids <- collated_seqone_csv_data_mod |> 
  inner_join(seqone_sample_info,
             by = "labno") |> 
  relocate(labno, firstname, surname, nhsno, dob)

stopifnot(anyNA(seqone_csv_data_with_ids$firstname) == FALSE)

# Remove validation data --------------------------------------------------

DOC6192_validation_worksheets <- c("WS133557", "WS134687", "WS134928", 
                                   "WS135001", "WS135498")

DOC6255_validation_worksheets <- c("WS136827", "WS138201", "WS138439", 
                                   "WS138627")

# 2 validation samples were on WS138061. The other samples on this
# worksheet were live clinical samples
DOC6255_validation_samples <- c("WS138061_23047082", "WS138061_23053359")

DOC6588_validation_worksheets <- c("WS147582", "WS149085", "WS149086")

validation_worksheets <- c(DOC6192_validation_worksheets,
                           DOC6255_validation_worksheets,
                           DOC6588_validation_worksheets)

seqone_data_val_removed <- seqone_csv_data_with_ids |> 
  filter(!worksheet %in% validation_worksheets) |> 
  filter(!sample %in% DOC6255_validation_samples)

# Remove reference control data -------------------------------------------

seqone_data_controls_removed <- seqone_data_val_removed |> 
  filter(!surname %in% c("Seraseq", "GenQA"))

# Remove repeat tests of same sample --------------------------------------

seqone_data_dups_removed <- seqone_data_controls_removed |> 
  filter(!duplicated(labno)) 

# Remove recent results ---------------------------------------------------

# Remove any results after 21st April 2025 as these have not been 
# reported yet.

seqone_data_recent_removed <- seqone_data_dups_removed |> 
  filter(date < "2025-04-22  UTC")

# Remove multiple tests from same patient ---------------------------------

seqone_data_unique_patient_results <- seqone_data_recent_removed |> 
  filter(!is.na(nhsno)) |> 
  # Some patient have multiple samples. To select samples with conclusive 
  # results (positive or negative statuses), arrange by NHS number and 
  # status and then
  # remove NHS number duplicates
  arrange(nhsno, status) |> 
  filter(!duplicated(nhsno))

seqone_data_patient_repeat_results <- seqone_data_recent_removed |> 
  filter(!is.na(nhsno)) |> 
  filter(duplicated(nhsno, fromLast = TRUE) |
           duplicated(nhsno, fromLast = FALSE)) |> 
  arrange(nhsno, status, labno) |> 
  select(labno, firstname, surname, nhsno, score, status)

# I've checked and only one patient within the repeated results has only
# non-conclusive results
seqone_patients_with_only_non_conclusive_results <- seqone_data_unique_patient_results |> 
  filter(labno %in% seqone_data_patient_repeat_results$labno &
           status == "Non-conclusive")

stopifnot(seqone_patients_with_only_non_conclusive_results$labno == "24009804")

# For samples without NHS numbers, I manually checked them to make sure that 
# none of these samples were from the same patient based on firstname and
# surname
seqone_data_without_nhsno <- seqone_data_dups_removed |> 
  filter(is.na(nhsno))

esmo_data <- rbind(seqone_data_without_nhsno, 
                   seqone_data_unique_patient_results) 

stopifnot(anyDuplicated(esmo_data$nhsno, incomparables = c(NA)) == 0)

stopifnot(anyDuplicated(esmo_data$labno) == 0)

# Remove patient identifiers ----------------------------------------------

esmo_data_ids_removed <- esmo_data |> 
  select(labno, worksheet, sample, date, 
         somahrd_version, LGA, LPC, score, status,
         ccne1_cn, rad51b_cn, coverage, pct_mapped_reads,
         gi_confidence,
         low_tumor_fraction)

# Export collated data ----------------------------------------------------

write_csv(esmo_data_ids_removed,
          file = paste0(esmo_folderpath,
                        "collated/esmo_data_ids_removed.csv"))
