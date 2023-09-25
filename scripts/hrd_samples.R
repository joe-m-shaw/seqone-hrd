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

source("scripts/hrd_filepaths.R")
source("scripts/dlms_connection.R")

#############################
# DNA Volumes
#############################

hrd_sample_volumes <- read_csv(paste0(hrd_data_path, 
                                      "hrd_sample_volumes.csv")) %>%
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
dlms_204 <- read_csv(paste0(hrd_data_path, "DDBK_Samples_204.csv")) %>%
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

hrd_sample_ncc <- rbind(dlms_204 %>%
                          select(dlms_dna_number, ncc),
                        sample_info_additional_samples %>%
                          select(dlms_dna_number, ncc)) %>%
  filter(!base::duplicated(dlms_dna_number)) %>%
  mutate(dlms_dna_number = as.numeric(dlms_dna_number))

samples_needing_ncc <- hrd_sample_ncc %>%
  filter(is.na(ncc))

# NCC information checked by Katie Sadler
ncc_ks <- read_csv(paste0(hrd_data_path, 
                          "annotated_hrd_sample_list_KS.csv")) %>%
  janitor::clean_names() %>%
  dplyr::rename(dlms_dna_number = specimen_number)

ncc_ks_clean <- ncc_ks %>%
  filter(dlms_dna_number %in% samples_needing_ncc$dlms_dna_number) %>%
  select(dlms_dna_number, checked_ncc) %>%
  dplyr::rename(ncc = checked_ncc)

hrd_sample_ncc <- rbind(hrd_sample_ncc, ncc_ks_clean) %>%
  filter(!is.na(ncc))

# Extra DLMS export for code 215 (DNA store) - includes information
# for biobank and seraseq controls
dlms_215 <- read_csv(paste0(hrd_data_path, "DDBK_Samples_215.csv")) %>%
  janitor::clean_names() %>%
  dplyr::rename(dlms_dna_number = labno) 


dlms_joined <- rbind(dlms_204 %>%
                       select(-c(ncc_single,
                                 ncc_range,
                                 ncc)), dlms_215)

#############################
# DNA concentrations
#############################

# Sample concentrations sent by Rebecca Hall, 
# collated from sequencing spreadsheets
hrd_sample_concentrations <- read_csv(paste0(hrd_data_path, 
                                         "hrd_sample_concentrations.csv")) %>%
  janitor::clean_names() %>%
  dplyr::rename(dlms_dna_number = specimen_number) %>%
  select(dlms_dna_number, qubit_ng_u_l) %>%
  mutate(qubit_ng_u_l = as.numeric(qubit_ng_u_l))

#############################
# Previous SeqOne Runs
#############################

seqone_run1 <- read_csv(paste0(hrd_data_path, "seqone_run1.csv")) %>%
  filter(!base::duplicated(specimen_number))

#############################
# Pathology Block IDs
#############################

pathology_block_ids <- read_csv(paste0(hrd_data_path, "pathology_block_ids.csv"))

##################################################
# Join Sample Information
##################################################

colnames(dlms_joined)

collated_hrd_sample_info <- hrd_sample_volumes %>%
  left_join(hrd_sample_gi_scores, by = "dlms_dna_number") %>%
  left_join(hrd_sample_concentrations, by = "dlms_dna_number") %>%
  left_join(hrd_sample_ncc, by = "dlms_dna_number") %>%
  left_join(dlms_joined %>%
              select(dlms_dna_number, i_gene_r_no, firstname, surname,
                     nhsno), 
            by = "dlms_dna_number") %>%
  mutate(seqone_run1 = ifelse(dlms_dna_number %in% seqone_run1$specimen_number, 
         "Yes", "No"),
         gi_score = as.numeric(gi_score),
         
         myriad_hrd_result = case_when(
              gi_score >= 42 ~"Positive",
              gi_score < 42 ~"Negative")) %>%
  dplyr::rename(nhs_number = nhsno)

##################################################
# Selecting samples for run 2 and run 3
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

run_three_samples_non_repeat <- c(23016922, 23016919, 23016917, 23016815,
                                   23016526, 23016524, 23016520, 23016518,
                                   23016516, 21022554, 21014781, 21006236,
                                   21002698, 21000693, 20119814,
                                   20119049, 20113426, 20110258, 23017428,
                                   23025928, 23033359)

run_three_samples <- c(run_three_samples_non_repeat, repeat_samples,
                       # Accidentally repeated this sample on run 2, so will
                       # repeat on run3 
                       21013520)

annotated_hrd_sample_list <- collated_hrd_sample_info %>%
  mutate(
    seqone_run2 = ifelse(dlms_dna_number %in% run_two_samples, "Yes", "No"),
    
    seqone_run3 = ifelse(dlms_dna_number %in% run_three_samples, "Yes", "No"),
    
    reason = case_when(
      dlms_dna_number %in% low_quality_samples ~"New sample: poor quality",
      dlms_dna_number %in% repeat_samples ~"Repeat from run one",
      dlms_dna_number %in% new_samples ~"New sample: better quality")) %>%
  select(dlms_dna_number, i_gene_r_no, firstname, surname, nhs_number, 
         volume_ul, qubit_ng_u_l, 
         ncc, gi_score, 
         myriad_hrd_result, seqone_run1, seqone_run2, seqone_run3) %>%
  arrange(desc(seqone_run1), desc(seqone_run2), desc(seqone_run3)) %>%
  mutate()


test2 <- get_sample_data(annotated_hrd_sample_list$dlms_dna_number)

#############################
# Spreadsheet from Katie Sadler
#############################

sample_brca_spreadsheet <- read_excel(paste0(hrd_data_path, 
           "HRD TBRCA data collection Manchester_NEW_from Oct2022_2023.xlsx"),
    skip  = 1) %>%
  janitor::clean_names() %>%
  dplyr::rename(i_gene_r_no = lab_number_id) %>%
  mutate(gis_numeric = as.numeric(gis_score_numerical_value_or_fail_not_tested))

samples_to_add <- dlms_204 %>%
  left_join(sample_brca_spreadsheet, by = "i_gene_r_no") %>%
  filter(!is.na(gis_numeric) & 
           !dlms_dna_number %in% annotated_hrd_sample_list$dlms_dna_number)

new_samples_for_volume_check <- rbind(samples_to_add %>%
  select(dlms_dna_number),
  new_samples_katie %>%
    select(dlms_dna_number))

##################################################
# Spread of Genomic Instability Scores
##################################################

samples_for_run4_from_run_3 <- c(20112141, 
                                  23016815, 
                                  23016917, 
                                  20104105,
                                  23016919,
                                  21022554,
                                  21014781)

samples_for_run4_from_run_2 <- c(20112141,
                                 21012359,
                                 21013520,
                                 21006928,
                                 21009404,
                                 21008398,
                                 21009934)

samples_to_test <- annotated_hrd_sample_list %>%
  filter(dlms_dna_number %in% samples_for_run4_from_run_2) 

ggplot(samples_to_test, aes(x = reorder(dlms_dna_number, gi_score), 
                            y = gi_score,
                            colour = myriad_hrd_result)) +
  geom_point(size = 2) +
  ylim(0, 100) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90)) +
  labs(x = "", y = "Myriad Genomic Instability Score",
       title = "Variation in Myriad Genomic Instability Scores") +
  geom_hline(yintercept = 42, linetype = "dashed")


write.csv(samples_to_test %>%
            select(dlms_dna_number, i_gene_r_no,
                   firstname, surname),
          paste0(hrd_project_path, "outputs/7_samples_fo-next_run.csv"),
          row.names = FALSE)


WS134687 <- annotated_hrd_sample_list %>%
  filter(seqone_run2 == "Yes") %>%
  select(dlms_dna_number, firstname, surname,
         ncc)

write.csv(WS134687,paste0(hrd_project_path, "outputs/WS134687_ncc.csv"),
          row.names = FALSE)

##################################################
# Sample Concentrations
##################################################

# Check how many samples would be excluded by setting DNA input limits
# at different thresholds.

dlms_204 %>%
  filter(concentration < 500) %>%
  ggplot(aes(x = disease, y = concentration)) +
  geom_boxplot() +
  theme_bw() +
  theme(panel.grid = element_blank())

# 50ng input at 15ul
upper_threshold <- 50/15
# 15ng input at 15ul
lower_threshold <- 15/15

check_thresholds <- dlms_204 %>%
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
# Export
##################################################

write.csv(annotated_hrd_sample_list, 
          file = paste0(hrd_project_path, "outputs/annotated_hrd_sample_list_", 
                        format(Sys.time(), "%Y_%m_%d_%H_%M_%S"),
                        ".csv"),
          row.names = FALSE)

##################################################
# Samples with BRCA variants for PanSolid Run
##################################################

samples_for_pansolid <- c(23016917, 23016516, 23016815,
                          21001739, 21006858)


samples_for_pansolid_details <- annotated_hrd_sample_list %>%
  filter(dlms_dna_number %in% samples_for_pansolid) %>%
  select(dlms_dna_number, i_gene_r_no, firstname, surname)

write.csv(samples_for_pansolid_details, 
          paste0(hrd_project_path, "outputs/samples_for_pansolid_details.csv"),
          row.names = FALSE)

##################################################