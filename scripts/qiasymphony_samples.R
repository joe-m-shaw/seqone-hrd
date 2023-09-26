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

##################################################
# Input
##################################################

input_samples <- read_csv(paste0(hrd_data_path, "louise_samples.csv"))

##################################################
# Extract data from DLMS tables
##################################################

sample_info <- get_sample_data(input_samples$lab_no) %>%
  select(labno, disease, disease_2, disease_3, disease_4, pathno,
         comments, concentration, nanodrop_ratio, nanodrop230ratio)
  
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

# Currently doesn't work
# gene_info <- sqlQuery(channel = moldb_connection,
                      #query = "SELECT * FROM MolecularDB.dbo.ResultsAccess WHERE LABNO IN (23039903)")

##################################################
# Join and export
##################################################

output <- sample_info %>%
  left_join(extraction_batch_info, by = "labno") %>%
  left_join(extraction_batch_dates, by = "extraction_batch_id") %>%
  arrange(labno) %>%
  filter(extraction_method_fk == 25)
  
write.csv(output, paste0(hrd_project_path, "outputs/qia_symphony_extraction_batch_info.csv"),
          row.names = FALSE)

##################################################
