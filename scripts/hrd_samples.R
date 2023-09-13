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

# Sample information exported from DLMS with disease code 204
dlms_info <- read_csv(paste0(homepath, "data/DDBK_Samples_204.csv")) %>%
  janitor::clean_names() %>%
  dplyr::rename(specimen_number = labno)

##################################################
# Join Sample Information
##################################################

collated_sample_info <- sample_volumes %>%
  dplyr::rename(volume_ul_12_09_2023 = volume_ul) %>%
  mutate(seqone_run1 = ifelse(specimen_number %in% seqone_run1$specimen_number, 
         "Yes", "No")) %>%
  left_join(dlms_info, by = "specimen_number")

##################################################
# Extract neoplastic cell content using regex
##################################################

# Can be NCC or TCC
comment_regex_single <- ".+(\\W{1}\\d{2}%).+"

comment_regex_range <- ".+(\\d{2}\\W{1}\\d{2}%).+"

check <- collated_sample_info %>%
  mutate(ncc_single = sub(x = comments,
                   pattern = comment_regex_single,
                   replacement = "\\1"),
         ncc_range = sub(x = comments,
                          pattern = comment_regex_range,
                          replacement = "\\1"),
         ncc_final = ifelse(length(ncc_range) == 6, ncc_range, ncc_single)) %>%
  select(specimen_number, comments, ncc_single, ncc_range, ncc_final)

##################################################
# Export
##################################################

#write.csv(collated_sample_info, 
          #file = paste0(homepath, "outputs/collated_hrd_information.csv"),
          #row.names = FALSE)

##################################################