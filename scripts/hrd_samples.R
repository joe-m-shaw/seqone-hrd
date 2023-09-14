################################################################################
# SeqOne Samples
# joseph.shaw3@nhs.net
################################################################################

##################################################
# Packages
##################################################

library(tidyverse)
library(readxl)
library(janitor)

##################################################
# Load Sample Information
##################################################

homepath <- "~/homologous_recombination_deficiency/"

sample_volumes <- read_csv(paste0(homepath, "data/hrd_sample_volumes.csv")) %>%
  janitor::clean_names()

seqone_run1 <- read_csv(paste0(homepath, "data/seqone_run1.csv")) %>%
  filter(!base::duplicated(specimen_number))

# Sample spreadsheet sent by Eleanor Baker
sample_info <- read_excel(paste0(homepath, "data/HRD Trial samples - QCs added(AutoRecovered).xlsx"),
                          sheet = "Sheet1") %>%
  janitor::clean_names() %>%
  mutate(specimen_number = as.numeric(lab_no))

# Can be NCC or TCC
comment_regex_single <- ".+(\\W{1}\\d{2}%).+"

comment_regex_range <- ".+(\\d{2}\\W{1}\\d{2}%).+"

# Sample information exported from DLMS with disease code 204
dlms_info <- read_csv(paste0(homepath, "data/DDBK_Samples_204.csv")) %>%
  janitor::clean_names() %>%
  dplyr::rename(specimen_number = labno) %>%
  # Extract neoplastic cell content using regex
  mutate(ncc_single = sub(x = comments,
                          pattern = comment_regex_single,
                          replacement = "\\1"),
         ncc_range = sub(x = comments,
                         pattern = comment_regex_range,
                         replacement = "\\1"),
         ncc_final = ifelse(str_length(ncc_single) == 4 & ncc_range > 6,
                            ncc_single, ifelse(str_length(ncc_range) == 6, ncc_range, NA)))

# Sample concentrations sent by Rebecca Hall
sample_concentrations <- read_csv(paste0(homepath, "data/hrd_sample_concentrations.csv")) %>%
  janitor::clean_names() %>%
  select(specimen_number, qubit_ng_u_l) %>%
  mutate(qubit_ng_u_l = as.numeric(qubit_ng_u_l))

# Sample spreadsheet sent by Katie Sadler
sample_brca_spreadsheet <- read_excel(paste0(homepath, "data/HRD TBRCA data collection Manchester_NEW_from Oct2022_2023.xlsx"),
                                      skip  = 1) %>%
  janitor::clean_names() %>%
  dplyr::rename(i_gene_r_no = lab_number_id)

colnames(sample_brca_spreadsheet)

##################################################
# Join Sample Information
##################################################

collated_sample_info <- sample_volumes %>%
  dplyr::rename(volume_ul_12_09_2023 = volume_ul) %>%
  mutate(seqone_run1 = ifelse(specimen_number %in% seqone_run1$specimen_number, 
         "Yes", "No")) %>%
  left_join(dlms_info, by = "specimen_number") %>%
  left_join(sample_concentrations, by = "specimen_number") %>%
  select(-concentration) %>%
  mutate(tests_remaining = round((volume_ul_12_09_2023 * qubit_ng_u_l) / 50)) %>%
  left_join(sample_info %>%
              select(specimen_number, gis_score), by = "specimen_number") %>%
  select(i_gene_r_no, specimen_number, volume_ul_12_09_2023, qubit_ng_u_l, tests_remaining,
         seqone_run1, gis_score, ncc_final)
        
##################################################
# Export
##################################################

#write.csv(collated_sample_info, 
          #file = paste0(homepath, "outputs/collated_hrd_information.csv"),
          #row.names = FALSE)

##################################################
# Selecting samples for run 2
##################################################

# Quarter of samples focussed on reproducibility (8)
# Third on pushing limits – less than 50ng (10)
# Rest: range of GIS scores, with better input (13)

##################################################