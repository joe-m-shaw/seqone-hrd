# Select samples for manual library prep validation

library(tidyverse)

# Cohort 1: Cobas-extracted DNA from the initial validation repeat study

seqone_validation_data <- read_csv(file = paste0(config::get("data_filepath"), 
                                               "validation/DOC6192_seqone/",
                                               "outputs/",
                                               "2023_11_10_10_55_08_compare_results.csv"),
                                   col_types = list(
                                     "dlms_dna_number" = col_character(),
                                     "version" = col_character()
                                   )) 

repeat_samples <- c("20112141", "20103853", "21003549",
                    "20127786", "21011999", "21012359",
                    "20103853", 
                    # Low DNA: false negative in the validation
                    "23016526")

cohort1 <- seqone_validation_data |> 
  filter(dlms_dna_number %in% repeat_samples &
           version == "1.2") |> 
  select(dlms_dna_number, firstname, surname, lga, lpc, seqone_hrd_score, 
         seqone_hrd_status, myriad_r_number, myriad_gi_score, myriad_hrd_status,
         path_block_manual_check, qubit_dna_ul) |> 
  filter(!duplicated(dlms_dna_number)) |> 
  mutate(cohort = "Cohort 1")

# Cohort 2: samples with borderline SeqOne HRD scores

seqone_live_csv_data <- read_csv(file = paste0(config::get("data_filepath"), 
                                               "live_service/service/",
                                               "collated_data/",
                                               "collated_live_csv_data.csv"),
                                 col_types = "cccdddclldddddddcc") |> 
  janitor::clean_names()

cohort2 <- seqone_live_csv_data |> 
  mutate(labno = str_extract(string = sample, pattern = "\\d{8}")) |> 
  filter(labno %in% c("24055693", "24021526", "24056470", 
                      "24057137", "24035637", "24015562",
                      # Inconclusive sample
                      "24063655")) |> 
  mutate(cohort = "Cohort 2") |> 
  select(labno, lga, lpc, score, status)


write_csv(cohort1, "cohort1.csv")

write_csv(cohort2, "cohort2.csv")
