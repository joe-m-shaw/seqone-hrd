################################################################################
# QiaSymphony Validation Samples
# joseph.shaw3@nhs.net
################################################################################

rm(list=ls())

##################################################
# Source scripts
##################################################

source("scripts/hrd_filepaths.R")
source("scripts/dlms_connection.R")

library(tidyverse)

##################################################
# Input
##################################################

input_samples <- read_csv(paste0(hrd_data_path, "louise_samples.csv"))

##################################################
# Extract data from DLMS tables
##################################################

sample_info <- get_sample_data(input_samples$lab_no) %>%
  select(labno, dob, firstname, surname, i_gene_r_no, i_gene_s_no, 
         nhsno, disease, pathno, comments, concentration, nanodrop_ratio, 
         nanodrop230ratio) %>%
  mutate(dob_mod = sub(pattern = "^(\\d{4}-\\d{2}-\\d{2}).+",
                       replacement = "\\1",
                       x =dob))
  
extraction_batch_info <- sqlQuery(channel = moldb_connection,
                          query = paste0("SELECT * FROM MolecularDB.dbo.MOL_Extractions WHERE LABNO IN (",
                                         paste(input_samples$lab_no, collapse = ", "),
                                         ")")) %>%
  arrange(LabNo) %>%
  janitor::clean_names() %>%
  dplyr::rename(labno = lab_no,
                extraction_batch_id = extraction_batch_fk)

extraction_batches <- unique(extraction_batch_info$extraction_batch_id)

extraction_batch_dates <- sqlQuery(channel = moldb_connection,
         query = paste0("SELECT * FROM MolecularDB.dbo.MOL_ExtractionBatches WHERE ExtractionBatchId IN (",
                        paste(extraction_batches, collapse = ", "),
                        ")")) %>%
  janitor::clean_names()

# Workaround - selecting by LABNO doesn't work due to lack of data consistency
# in LABNO column. Instead, pull out everything via plate position, then filter.
input_numbers <- seq(1, 96, 1)

results_access <- sqlQuery(channel = moldb_connection,
                      query = paste0("SELECT LABNO, DISCODE, Exon, Genotype FROM MolecularDB.dbo.ResultsAccess WHERE POSITION IN (", 
                 paste(input_numbers, collapse = ", "),
                 ")")) %>%
  janitor::clean_names()

exon_table <- results_access %>%
  filter(labno %in% input_samples$lab_no) %>%
  # Change data type to allow join in next step
  mutate(labno = as.numeric(labno))

##################################################
# Join
##################################################

output <- sample_info %>%
  left_join(extraction_batch_info, by = "labno") %>%
  left_join(extraction_batch_dates, by = "extraction_batch_id") %>%
  left_join(exon_table, by = "labno", relationship = "many-to-many") %>%
  arrange(labno) %>%
  select(-c(status, macro, version.x, version.y, 
            extraction_id, locked, dob))

##################################################
# Check
##################################################

stopifnot(setdiff(input_samples$lab_no, output$labno) == 0)

##################################################
# Export
##################################################

louise_folder <- "S:/central shared/Genetics/Mol_Shared/Louise.Kung/"

write.csv(output, paste0(louise_folder, "qia_symphony_extraction_batch_info",
                         format(Sys.time(), "%Y_%m_%d_%H_%M_%S"),
                         ".csv"),
          row.names = FALSE)

##################################################