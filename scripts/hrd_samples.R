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

hrd_project_path <- "~/homologous_recombination_deficiency/"

hrd_data_path <- "~/homologous_recombination_deficiency/data/"

#############################
# DNA Volumes
#############################

hrd_sample_volumes <- read_csv(paste0(hrd_data_path, "hrd_sample_volumes.csv")) %>%
  janitor::clean_names() %>%
  dplyr::rename(dlms_dna_number = specimen_number)

#############################
# Genomic Instability Scores
#############################

# Sample spreadsheet sent by Eleanor Baker
sample_info <- read_excel(paste0(hrd_data_path, 
                         "HRD Trial samples - QCs added(AutoRecovered).xlsx"),
                  sheet = "Sheet1") %>%
  janitor::clean_names() %>%
  mutate(dlms_dna_number = as.numeric(lab_no)) %>%
  dplyr::rename(gi_score = gis_score)

sample_info_additional_samples <- read_excel(paste0(hrd_data_path, 
                          "HRD Trial samples - QCs added(AutoRecovered).xlsx"),
                   sheet = "Additional samples") %>%
  janitor::clean_names() %>%
  mutate(dlms_dna_number = as.numeric(lab_no)) %>%
  filter(!is.na(lab_no) & lab_no != "Not found" & !is.na(gis)) %>%
  dplyr::rename(gi_score = gis)

# Samples from .txt document sent by Katie Sadler
new_samples_katie <- read_csv(paste0(hrd_data_path, 
                                     "samples_from_katie_sadler.csv")) %>%
  mutate(gis_numeric = as.numeric(gsub("\\D", "", gis_score))) %>%
  filter(!is.na(gis_numeric)) %>%
  dplyr::rename(gi_score = gis_numeric,
                dlms_dna_number = specimen_number)

hrd_sample_gi_scores <- rbind(sample_info %>%
                              select(dlms_dna_number, gi_score),
                             sample_info_additional_samples %>%
                               select(dlms_dna_number, gi_score),
                             new_samples_katie %>%
                               select(dlms_dna_number, gi_score))

#############################
# Neoplastic Cell Content
#############################

# Tumour cell content can be input as "NCC" or "TCC"
# and as a single value with an operator (example ">30%")
# or as a range (example: "10-20%")
comment_regex_single <- ".+(\\W{1}\\d{2}%).+"

comment_regex_range <- ".+(\\d{2}\\W{1}\\d{2}%).+"

# Sample information exported from DLMS with disease code 204
dlms_info <- read_csv(paste0(hrd_data_path, "DDBK_Samples_204.csv")) %>%
  janitor::clean_names() %>%
  dplyr::rename(dlms_dna_number = labno) %>%
  # Extract neoplastic cell content using regex
  mutate(ncc_single = sub(x = comments,
                          pattern = comment_regex_single,
                          replacement = "\\1"),
         ncc_range = sub(x = comments,
                         pattern = comment_regex_range,
                         replacement = "\\1"),
         ncc = ifelse(str_length(ncc_single) == 4 & ncc_range > 6,
                            ncc_single, ifelse(str_length(ncc_range) == 6, 
                                               ncc_range, NA))) %>%
  filter(!base::duplicated(i_gene_r_no))

hrd_sample_ncc <- rbind(dlms_info %>%
                          select(dlms_dna_number, ncc),
                        sample_info_additional_samples %>%
                          select(dlms_dna_number, ncc)) %>%
  filter(!base::duplicated(dlms_dna_number)) %>%
  mutate(dlms_dna_number = as.numeric(dlms_dna_number))

#############################
# DNA concentrations
#############################

# Sample concentrations sent by Rebecca Hall, collated from sequencing spreadsheets
hrd_sample_concentrations <- read_csv(paste0(hrd_data_path, 
                                         "hrd_sample_concentrations.csv")) %>%
  janitor::clean_names() %>%
  dplyr::rename(dlms_dna_number = specimen_number) %>%
  select(dlms_dna_number, qubit_ng_u_l) %>%
  mutate(qubit_ng_u_l = as.numeric(qubit_ng_u_l))

#############################
# Spreadsheet from Katie Sadler
#############################

sample_brca_spreadsheet <- read_excel(paste0(homepath, 
        "data/HRD TBRCA data collection Manchester_NEW_from Oct2022_2023.xlsx"),
                                      skip  = 1) %>%
  janitor::clean_names() %>%
  dplyr::rename(i_gene_r_no = lab_number_id) %>%
  mutate(gis_numeric = as.numeric(gis_score_numerical_value_or_fail_not_tested))

samples_to_add <- dlms_info %>%
  left_join(sample_brca_spreadsheet, by = "i_gene_r_no") %>%
  filter(!is.na(gis_numeric) & !specimen_number %in% annotated_hrd_sample_list$specimen_number)

#############################
# Previous SeqOne Runs
#############################

seqone_run1 <- read_csv(paste0(hrd_data_path, "seqone_run1.csv")) %>%
  filter(!base::duplicated(specimen_number))

##################################################
# Join Sample Information
##################################################

collated_hrd_sample_info <- hrd_sample_volumes %>%
  left_join(hrd_sample_gi_scores, by = "dlms_dna_number") %>%
  left_join(hrd_sample_concentrations, by = "dlms_dna_number") %>%
  left_join(hrd_sample_ncc, by = "dlms_dna_number") %>%
  mutate(seqone_run1 = ifelse(dlms_dna_number %in% seqone_run1$specimen_number, 
         "Yes", "No")) %>%
  mutate(myriad_hrd_result = case_when(
    gi_score >= 42 ~"Positive",
    gi_score < 42 ~"Negative"))

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
  filter(!is.na(gi_score)) %>%
  filter(qubit_ng_u_l < 3.3)

random_samples <- collated_hrd_sample_info %>%
  filter(seqone_run1 == "No") %>%
  filter(!is.na(gi_score)) %>%
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
  21013520, 21012001, 21008471)

repeat_samples <- c(
  # Seraseq control samples
  23031639, 23032088, 23032086, 
  # Breast cancer patients
  20103853, 20104105, 20112141, 20127786, 21012359, 21003549)

run_two_samples <- c(low_quality_samples, new_samples, repeat_samples)

annotated_hrd_sample_list <- collated_hrd_sample_info %>%
  mutate(
    seqone_run2 = ifelse(dlms_dna_number %in% run_two_samples, "Yes", "No"),
    
    reason = case_when(
      dlms_dna_number %in% low_quality_samples ~"New sample: poor quality",
      dlms_dna_number %in% repeat_samples ~"Repeat from run one",
      dlms_dna_number %in% new_samples ~"New sample: better quality")) %>%
  arrange(reason)

# Draw a plot
annotated_hrd_sample_list %>%
  filter(specimen_number %in% low_quality_samples |
           specimen_number %in% new_samples 
         #| specimen_number %in% repeat_samples
         ) %>%
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

write.csv(annotated_hrd_sample_list %>%
            filter(seqone_run2 == "Yes"), 
          file = paste0(homepath, "outputs/seqone_run2_samples_", 
                                         format(Sys.time(), "%Y_%m_%d_%H_%M_%S"),
                                         ".csv"),
          row.names = FALSE)

##################################################
# Selecting samples for run 3
##################################################

# We need 22 samples
# 4 from Katie's spreadsheet
# 12 from the annotated spreadsheet
# 6

##################################################
# Spread of Genomic Instability Scores
##################################################

samples_to_test <- annotated_hrd_sample_list %>%
  filter((seqone_run1 == "Yes" | seqone_run2 == "Yes") &
           !is.na(gis_score)) 

ggplot(samples_to_test, aes(x = reorder(specimen_number, gis_score), y = gis_score)) +
  geom_point(size = 2) +
  ylim(0, 100) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90)) +
  labs(x = "", y = "Myriad Genomic Instability Score")

check_gis <- annotated_hrd_sample_list %>%
  filter((seqone_run1 == "Yes" | seqone_run2 == "Yes") &
           !is.na(gis_score))  %>%
  select(specimen_number, gis_score)

max(samples_to_test$gis_score)
min(samples_to_test$gis_score)

##################################################
# Sample Concentrations
##################################################

# Check how many samples would be excluded by setting DNA input limits
# at different thresholds.

dlms_info %>%
  filter(concentration < 500) %>%
  ggplot(aes(x = disease, y = concentration)) +
  geom_boxplot() +
  theme_bw() +
  theme(panel.grid = element_blank())

# 50ng input at 15ul
upper_threshold <- 50/15
# 15ng input at 15ul
lower_threshold <- 15/15

check_thresholds <- dlms_info %>%
  filter(!is.na(concentration)) %>%
  mutate(
    pass_upper_threshold = case_when(
    concentration >= upper_threshold ~"Yes",
    concentration < upper_threshold ~"No"),
    
    pass_lower_threshold = case_when(
      concentration >= lower_threshold ~"Yes",
      concentration < lower_threshold ~"No"))

check_thresholds %>%
  group_by(pass_upper_threshold) %>%
  summarise(total = n()) %>%
  mutate(percent = (round(total / sum(total), 3))*100)

check_thresholds %>%
  group_by(pass_lower_threshold) %>%
  summarise(total = n()) %>%
  mutate(percent = (round(total / sum(total), 3))*100)

##################################################