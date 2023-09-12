################################################################################
# SeqOne Samples
# joseph.shaw3@nhs.net
################################################################################

rm(list=ls())

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

sample_info <- read_excel(paste0(homepath, "data/HRD Trial samples - QCs added(AutoRecovered).xlsx"),
                          sheet = "Sheet1") %>%
  janitor::clean_names() %>%
  mutate(specimen_number = as.numeric(lab_no))

##################################################
# Join Sample Information
##################################################

collated_sample_info <- sample_volumes %>%
  dplyr::rename(volume_ul_12_09_2023 = volume_ul) %>%
  mutate(seqone_run1 = ifelse(specimen_number %in% seqone_run1$specimen_number, 
         "Yes", "No")) %>%
  left_join(sample_info, by = "specimen_number") %>%
  select(-c(lab_no, sample_stratergy, date_reported, x12, x14, volume,
            conc, dqn))

write.csv(collated_sample_info, 
          file = paste0(homepath, "outputs/collated_hrd_information.csv"),
          row.names = FALSE)

##################################################