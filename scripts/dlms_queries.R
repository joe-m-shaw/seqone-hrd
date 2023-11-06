# DLMS queries

library("RODBC")

# Connection setup for PC38698
moldb_connection <- RODBC::odbcConnect(dsn = "moldb")

# SeqOne sample information ---------------------------------------------------------

seqone_dlms_info <- get_sample_data(collated_seqone_info$dlms_dna_number) |>
  rename(
    dlms_dna_number = labno,
    nhs_number = nhsno
  ) |>
  mutate(
    ncc_single = sub(
      x = comments,
      pattern = ".+(\\W{1}\\d{2}%).+",
      replacement = "\\1"
    ),
    ncc_range = sub(
      x = comments,
      pattern = ".+(\\d{2}\\W{1}\\d{2}%).+",
      replacement = "\\1"
    ),
    ncc = ifelse(str_length(ncc_single) == 4 & ncc_range > 6,
                 ncc_single, ifelse(str_length(ncc_range) == 6,
                                    ncc_range, NA
                 )
    )
  )

stopifnot(setdiff(
  seqone_dlms_info$dlms_dna_number,
  collated_seqone_info$dlms_dna_number
) == 0)

# Enter a fake NHS number for the Seraseq and Biobank controls, to allow
# Myriad scores to be joined later, and to keep Biobank controls within the
# dataset

seqone_dlms_info[seqone_dlms_info$dlms_dna_number == 23032086, "nhs_number"] <- 1
seqone_dlms_info[seqone_dlms_info$dlms_dna_number == 23032088, "nhs_number"] <- 2
seqone_dlms_info[seqone_dlms_info$dlms_dna_number == 23031639, "nhs_number"] <- 3
seqone_dlms_info[seqone_dlms_info$dlms_dna_number == 23033285, "nhs_number"] <- 4
seqone_dlms_info[seqone_dlms_info$dlms_dna_number == 23033279, "nhs_number"] <- 5
seqone_dlms_info[seqone_dlms_info$dlms_dna_number == 23033288, "nhs_number"] <- 6

write.csv(seqone_dlms_info, file = str_c(hrd_data_path, "seqone_dlms_info.csv"),
          row.names = FALSE)

# Script used to create csv for manual checking of pathology IDs

export_for_check <- seqone_dlms_info |>
  filter(dlms_dna_number %in% collated_seqone_info$dlms_dna_number) |>
  left_join(collated_myriad_info_mod |>
              select(
                nhs_number, myriad_pathology_block_pg1,
                myriad_pathology_block_pg2
              ), by = "nhs_number") |>
  select(
    dlms_dna_number, firstname, surname, nhs_number,
    pathno, myriad_pathology_block_pg1, myriad_pathology_block_pg2
  ) |>
  arrange(nhs_number)

# export_timestamp(hrd_data_path, export_for_check)

# tBRCA DLMS query ------------------------------------------------------------------

brca_query <- "SELECT * FROM MolecularDB.dbo.Samples WHERE DISEASE IN (204)"

tbrca_data <- sqlQuery(
  channel = moldb_connection,
  query = brca_query
) |>
  janitor::clean_names()

write.csv(tbrca_data, file = str_c(hrd_data_path, "tBRCA_dlms_info.csv"),
          row.names = FALSE)

# Sample extraction information -----------------------------------------------------

# Get tissues
seqone_tissue_types <- unique(seqone_dlms_info$tissue)

tissue_query <- sqlQuery(
  channel = moldb_connection,
  query = paste0(
    "SELECT * FROM MolecularDB.dbo.TissueTypes WHERE TissueTypeId IN (",
    paste(seqone_tissue_types, collapse = ", "),
    ")"
  )
) |>
  janitor::clean_names()

# Get extraction batch IDs
extraction_batch_id_table <- sqlQuery(
  channel = moldb_connection,
  query = paste0(
    "SELECT * FROM MolecularDB.dbo.MOL_Extractions WHERE LABNO IN (",
    paste(seqone_dlms_info$dlms_dna_number, collapse = ", "),
    ")"
  )
) |>
  arrange(LabNo) |>
  janitor::clean_names() |>
  dplyr::rename(
    dlms_dna_number = lab_no,
    extraction_batch_id = extraction_batch_fk
  )

extraction_batches <- unique(extraction_batch_id_table$extraction_batch_id)

# Get extraction batch dates
extraction_batch_date_table <- sqlQuery(
  channel = moldb_connection,
  query = paste0(
    "SELECT * FROM MolecularDB.dbo.MOL_ExtractionBatches WHERE ExtractionBatchId IN (",
    paste(extraction_batches, collapse = ", "),
    ")"
  )
) |>
  janitor::clean_names()

extraction_batch_info <- extraction_batch_id_table |>
  left_join(extraction_batch_date_table, by = "extraction_batch_id") |>
  # Get only extraction batch IDs used for Cobas extractions
  filter(extraction_method_fk == 25)

# 1 DNA number (23033279) is on 2 extraction batches

seqone_dlms_extractions <- seqone_dlms_info |>
  left_join(tissue_query |>
              dplyr::rename(tissue = tissue_type_id), by = "tissue") |>
  left_join(extraction_batch_info, by = "dlms_dna_number")

