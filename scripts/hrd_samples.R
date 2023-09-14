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
sample_info <- read_excel(paste0(homepath, 
                                 "data/HRD Trial samples - QCs added(AutoRecovered).xlsx"),
                          sheet = "Sheet1") %>%
  janitor::clean_names() %>%
  mutate(specimen_number = as.numeric(lab_no))

# Tumour cell content can be input as "NCC" or "TCC"
# and as a single value with an operator (example ">30%")
# or as a range (example: "10-20%")
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
                            ncc_single, ifelse(str_length(ncc_range) == 6, 
                                               ncc_range, NA))) %>%
  filter(!base::duplicated(i_gene_r_no))

# Sample concentrations sent by Rebecca Hall, collated from sequencing spreadsheets
sample_concentrations <- read_csv(paste0(homepath, 
                                         "data/hrd_sample_concentrations.csv")) %>%
  janitor::clean_names() %>%
  select(specimen_number, qubit_ng_u_l) %>%
  mutate(qubit_ng_u_l = as.numeric(qubit_ng_u_l))

# Sample spreadsheet sent by Katie Sadler
sample_brca_spreadsheet <- read_excel(paste0(homepath, 
        "data/HRD TBRCA data collection Manchester_NEW_from Oct2022_2023.xlsx"),
                                      skip  = 1) %>%
  janitor::clean_names() %>%
  dplyr::rename(i_gene_r_no = lab_number_id) %>%
  filter(!base::duplicated(i_gene_r_no))

##################################################
# Join Sample Information
##################################################

collated_hrd_sample_info <- sample_volumes %>%
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
         seqone_run1, gis_score, ncc_final) %>%
  mutate(myriad_hrd_result = case_when(
    gis_score >= 42 ~"Positive",
    gis_score < 42 ~"Negative"))

##################################################
# Selecting samples for run 2
##################################################

# Sample selection criteria:
# - Ovarian cancer patients
# - Have a result from Myriad HRD testing (Genomic Instability Score)
# - DNA from FFPE sample
# - Neoplastic cell content over 20%
# - Plenty of sample remaining for future testing

# Aim:
# 8 samples from run one, repeated for reproducibility
# 10 samples with less than 3.3ng/ul (less than 50ng input at 15ul), 
# to push the limits of the assay
# 13 samples with a range of GIS scores and better DNA concentrations

lower_quality_samples <- collated_hrd_sample_info %>%
  filter(seqone_run1 == "No") %>%
  filter(!is.na(gis_score)) %>%
  filter(qubit_ng_u_l < 3.3)

random_samples <- collated_hrd_sample_info %>%
  filter(seqone_run1 == "No") %>%
  filter(!is.na(gis_score)) %>%
  filter(qubit_ng_u_l > 3.3)

run_one_samples <- collated_hrd_sample_info %>%
  filter(seqone_run1 == "Yes")

low_quality_samples <- c(
  21009934, 21009404, 21006723, 21002912, 21001739, 20128167, 
  20126826)

new_samples <- c(
  # Negative GIS
  21016510, 21015168, 21008398, 21003078,
  20125849, 20112754, 21003752, 21016766,
  21035096, 21001595, 21006928, 21011999, 
  # Positive GIS
  21006858, 21012001, 21008471)

repeat_samples <- c(
  # Biobank samples
  23033288, 23033285, 23033279,
  # Breast cancer patients
  20103853, 20104105, 20112141, 20127786, 21012359, 21003549)

run_two_samples <- c(low_quality_samples, new_samples, repeat_samples)

annotated_hrd_sample_list <- collated_hrd_sample_info %>%
  mutate(
    seqone_run2 = ifelse(specimen_number %in% run_two_samples, "Yes", "No"),
    
    reason = case_when(
      specimen_number %in% low_quality_samples ~"New sample: poor quality",
      specimen_number %in% repeat_samples ~"Repeat from run one",
      specimen_number %in% new_samples ~"New sample: better quality")) %>%
  arrange(reason) %>%
  select(-c(i_gene_r_no, tests_remaining)) %>%
  dplyr::rename(volume_ul = volume_ul_12_09_2023,
                ncc = ncc_final)

# Draw a plot
annotated_hrd_sample_list %>%
  filter(specimen_number %in% low_quality_samples |
           specimen_number %in% new_samples |
           specimen_number %in% repeat_samples) %>%
  select(specimen_number, qubit_ng_u_l, gis_score) %>%
  pivot_longer(cols = -specimen_number,
               names_to = "category",
               values_to = "value") %>%
  mutate(specimen_number = as.character(specimen_number)) %>%
  ggplot(aes(x = specimen_number, y = value)) +
  geom_point(size = 2) +
  facet_wrap(~category) +
  labs(x = "Specimens") +
  theme_bw() +
  theme(axis.text.x = element_blank())

##################################################
# Export
##################################################

write.csv(annotated_hrd_sample_list, 
          file = paste0(homepath, "outputs/annotated_hrd_sample_list.csv"),
          row.names = FALSE)

##################################################