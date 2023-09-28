################################################################################
# Homologous Recombination Deficiency Validation: Analysis
# joseph.shaw3@nhs.net
################################################################################

rm(list=ls())

##################################################
# Packages and Filepaths
##################################################

library("ggpubr")
library("readxl")
library("epiR")

##################################################
# Source scripts
##################################################

source("scripts/collate_seqone_pdfs.R")
source("scripts/collate_myriad_pdfs.R")
source("scripts/dlms_connection.R")

downsampled_samples <- grep(pattern = "downsampl", 
                            x = collated_seqone_info$sample_id, value = TRUE)

##################################################
# Collate and edit Myriad Information
##################################################

# pdf_text doesn't work on these reports as resolution is too low.
# I had to collate data manually in an Excel
low_res_myriad_results <- read_excel(path = paste0(hrd_data_path,
                                                   "myriad_reports_low_res/myriad_low_res_report_data.xlsx")) %>%
  mutate(nhs_number = as.numeric(gsub(pattern = "(\\D)", "", nhs_number))) 

seraseq_gi_scores <- read_excel(paste0(hrd_data_path, "seraseq_gi_scores.xlsx")) %>%
  mutate(myriad_r_number = "",
         myriad_patient_name = "",
         myriad_dob = "",
         myriad_pathology_block = "",
         myriad_brca_status = "") %>%
  select(myriad_r_number, myriad_patient_name, myriad_dob, nhs_number,
         myriad_pathology_block, myriad_gi_score, myriad_hrd_status, 
         myriad_brca_status)

collated_myriad_info_mod <- rbind(collated_myriad_info, low_res_myriad_results,
                                  seraseq_gi_scores)

##################################################
# Pathology Block ID Check
##################################################

# Extract data from DLMS via ODBC
seqone_dlms_info <- get_sample_data(collated_seqone_info$dlms_dna_number) %>%
  dplyr::rename(dlms_dna_number = labno,
                nhs_number = nhsno) 

# Enter a fake NHS number for the Seraseq controls, to allow Myriad scores
# to be joined later.

seqone_dlms_info[seqone_dlms_info$dlms_dna_number == 23032086, "nhs_number"] <- 1
seqone_dlms_info[seqone_dlms_info$dlms_dna_number == 23032088, "nhs_number"] <- 2
seqone_dlms_info[seqone_dlms_info$dlms_dna_number == 23031639, "nhs_number"] <- 3

export_for_check <- collated_seqone_info %>%
  filter(!base::duplicated(dlms_dna_number)) %>%
  select(dlms_dna_number) %>%
  left_join(seqone_dlms_info %>%
              select(dlms_dna_number, firstname, surname, pathno, nhs_number), 
            by = "dlms_dna_number") %>%
  left_join(collated_myriad_info_mod %>%
              select(nhs_number, myriad_pathology_block), by = "nhs_number") %>%
  select(dlms_dna_number, firstname, surname, nhs_number, 
         pathno, myriad_pathology_block) %>%
  arrange(nhs_number) %>%
  mutate(path_block_automated_check = ifelse(pathno == myriad_pathology_block,
                                             "match", "NOT match"),
         path_block_manual_check = "")

# Writing of pathology block IDs is inconsistent, so a manual check is
# required

export_timestamp(hrd_data_path, export_for_check)

path_block_check <- read_csv(paste0(hrd_data_path, "manual_path_block_check_edit.csv")) %>%
  select(dlms_dna_number, path_block_manual_check)

##################################################
# Compare SeqOne and Myriad Results
##################################################

# Add QC data for SeqOne samples

seqone_qc_data <- read_excel(paste0(hrd_data_path, "seqone_qc_metrics_2023_09_28.xlsx")) %>%
  janitor::clean_names() %>%
  dplyr::rename(shallow_sample_id = sample,
                read_length = read_len,
                insert_size = ins_size,
                million_reads = m_reads)

seqone_mod <- collated_seqone_info %>%
  filter(!date %in% c("September 1, 2023", "August 31, 2023",
                      "August 24, 2023", "August 25, 2023")) %>%
  left_join(seqone_dlms_info %>%
              select(dlms_dna_number, nhs_number, firstname, surname, i_gene_r_no,
                     pathno),
            by = "dlms_dna_number") %>%
  mutate(downsampled = ifelse(sample_id %in% downsampled_samples, "Yes", "No")) %>%
  left_join(seqone_qc_data, by = "shallow_sample_id")

biobank_names <- unique(grep("Biobank", seqone_mod$surname,
                             value = TRUE, ignore.case = TRUE))

seraseq_names <- unique(grep("Seraseq", seqone_mod$surname,
                             value = TRUE, ignore.case = TRUE))

compare_results <- seqone_mod %>%
  left_join(collated_myriad_info_mod, by = "nhs_number") %>%
  left_join(path_block_check, by = "dlms_dna_number") %>%
  filter(!is.na(myriad_gi_score)) %>%
  mutate(hrd_status_check = case_when(
    
    myriad_hrd_status == "POSITIVE" & seqone_hrd_status == "POSITIVE" ~"Seqone HRD status consistent with Myriad",
    myriad_hrd_status == "NEGATIVE" & seqone_hrd_status == "NEGATIVE" ~"Seqone HRD status consistent with Myriad",
    myriad_hrd_status == "NEGATIVE" & seqone_hrd_status == "POSITIVE" ~"Seqone HRD status NOT consistent with Myriad",
    myriad_hrd_status == "POSITIVE" & seqone_hrd_status == "NEGATIVE" ~"Seqone HRD status NOT consistent with Myriad",
    TRUE ~"other"),
    
    identity = ifelse(surname %in% c(biobank_names, seraseq_names), "Control", "Patient"),
    control_type = case_when(
      surname %in% biobank_names ~"Biobank control",
      surname %in% seraseq_names ~"Seraseq control",
      TRUE ~"patient")) %>%
  filter(downsampled == "No")

compare_results[compare_results$dlms_dna_number == 23032086, "path_block_manual_check"] <- "pathology blocks match"
compare_results[compare_results$dlms_dna_number == 23032088, "path_block_manual_check"] <- "pathology blocks match"
compare_results[compare_results$dlms_dna_number == 23031639, "path_block_manual_check"] <- "pathology blocks match"

##################################################
# HRD Status Comparison
##################################################

positive_colour <- "#FF3333"
negative_colour <- "#3300FF"

comparison_summary <- compare_results %>%
  filter(!is.na(myriad_gi_score)) %>%
  group_by(hrd_status_check) %>%
  summarise(total = n())

ggplot(comparison_summary, aes(x = hrd_status_check, y = total)) +
  geom_col() +
  geom_text(aes(label = total), vjust = -0.5) +
  theme_bw() +
  labs(y = "Total DNA Inputs", x = "", title = "Consistency of Myriad and Seqone Testing")

##################################################
# Impact of Pathology Blocks
##################################################

results_for_path_block_plot <- compare_results %>%
  filter(!is.na(path_block_manual_check))

ggplot(results_for_path_block_plot, aes(x = myriad_gi_score, y = seqone_hrd_score)) +
  geom_point(size = 3, alpha = 0.6, aes(shape = hrd_status_check)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(limits = c(0,100),
                     breaks = c(0, 25, 42, 50, 75, 100)) +
  ylim(0, 1) +
  labs(x = "Myriad Genome Instability Score",
       y = "SeqOne HRD Score",
       title = "Comparison of Myriad vs SeqOne HRD Testing For Patient Samples",
       subtitle = paste0("Data for ", nrow(results_for_path_block_plot), " DNA inputs"),
       caption = "Repeat data included.") +
  #geom_vline(xintercept = 42, linetype = "dashed") +
  #geom_hline(yintercept = 0.5, linetype = "dashed") + 
  ggpubr::stat_cor(method = "pearson", label.x = 50, label.y = 0.25) +
  facet_wrap(~path_block_manual_check) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom",
        legend.title = element_blank())

##################################################
# Repeat Testing
##################################################

seqone_mod %>%
  filter(base::duplicated(dlms_dna_number, fromLast = TRUE) |
           base::duplicated(dlms_dna_number, fromLast = FALSE)) %>%
  filter(downsampled == "No") %>%
  ggplot(aes(x = worksheet, y = read_length)) +
  geom_point(size = 3, pch = 21) +
  facet_wrap(~dlms_dna_number) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "SeqOne results for repeated samples",
       x = "Worksheet",
       y = "Read length") +
  geom_hline(yintercept = 0.50, linetype = "dashed")

##################################################
# Impact of read length
##################################################

seqone_mod %>%
  filter(downsampled == "No") %>%
  ggplot(aes(x = reorder(shallow_sample_id, read_length), y = read_length)) +
  geom_point(size = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "SeqOne Read Length",
       x = "",
       y = "Read length")

##################################################
# Comparison with Simplified Model
##################################################

simplified_model <- compare_results %>%
  filter(myriad_brca_status == "NEGATIVE") %>%
  mutate(model_approximation = case_when(
    lga >= 18 ~"POSITIVE",
    lga <= 14 ~"NEGATIVE",
    lga < 18 & lga > 14 & lpc > 10 ~"POSITIVE",
    lga < 18 & lga > 14 & lpc <= 10 ~"NEGATIVE")) %>%
  mutate(approximation_consistent = ifelse(model_approximation == seqone_hrd_status,
                                        "Yes",
                                        "No"))

simplified_model_summary <- simplified_model %>%
  group_by(approximation_consistent) %>%
  summarise(total = n())

ggplot(simplified_model, aes(x = lga,
                            y = lpc,
                            fill = seqone_hrd_status)) +
  scale_fill_manual(values = c(negative_colour,
                               positive_colour)) +
  geom_point(size = 3, alpha = 0.6, aes(shape = approximation_consistent)) +
  scale_shape_manual(values = c(24, 21)) +
  theme_bw() +
  labs(x = "Large Genomic Alterations", 
       y = "Loss of Parental Copy",
       title = "Comparison of Seqone model approximation to pipeline output",
       subtitle = paste0("Data for ", nrow(simplified_model), " BRCA-negative DNA inputs")) 

##################################################
# Seraseq Controls
##################################################

seraseq_control_data <- compare_results %>%
  filter(control_type == "Seraseq control") %>%
  mutate(firstname_factor = factor(firstname, levels = c(
    "FFPE HRD Negative",
    "Low-Positive FFPE HRD",
    "High-Positive FFPE HRD")))

ggplot(seraseq_control_data, aes(x = myriad_gi_score, y = seqone_hrd_score)) +
  geom_point(size = 3, aes(shape = seqone_hrd_status,
                           colour = worksheet)) +
  ylim(0, 1) +
  xlim(0, 100) +
  facet_wrap(~firstname_factor) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "Myriad GI Score", y = "SeqOne HRD Score",
       title = "Seraseq Controls: repeat SeqOne data")

##################################################
# Export tables
##################################################

inconsistent_summary <- compare_results %>%
  filter(hrd_status_check == "Seqone HRD status NOT consistent with Myriad") %>%
  select(sample_id, worksheet, myriad_patient_name, seqone_hrd_score, myriad_gi_score,
         seqone_hrd_status, myriad_hrd_status, pathno, myriad_pathology_block, 
         path_block_manual_check, lga, lpc, ccne1, rad51b, coverage, percent_mapping) %>%
  arrange(sample_id, worksheet)

export_timestamp(hrd_output_path, inconsistent_summary)

unusual_samples <- c(21013520, 23032088, 20127786)

unusual_repeat_sample_info <- seqone_mod %>%
  filter(dlms_dna_number %in% unusual_samples) %>%
  select(-c(sample_id, filename, date, pathno, downsampled, user)) %>%
  arrange(dlms_dna_number)

export_timestamp(hrd_output_path, unusual_repeat_sample_info)

# Entire table of results
write.csv(compare_results, "S:/central shared/Genetics/Mol_Shared/Development.Team/SeqOne Homologous Recombination Deficiency Validation/2023_09_27 HRD validation update/result_comparison.csv",
          row.names = FALSE)

##################################################
# Save plots
##################################################

ggsave(plot = comparison_plot, 
       filename = paste0("comparison_plot_",
                         format(Sys.time(), "%Y%m%d_%H%M%S"),
                         ".jpg"),
       path = paste0(hrd_project_path, "plots/"), 
       device='jpeg',
       units = "cm",
       width = 15,
       height = 15)

##################################################
# Mid Output Run
##################################################

mid_output_run <- seqone_mod %>%
  filter(worksheet == "WS134928")

seqone_mod %>%
  filter(dlms_dna_number %in% mid_output_run$dlms_dna_number) %>%
  ggplot(aes(x = worksheet, y = lga)) +
  geom_point() +
  facet_wrap(~dlms_dna_number) 

##################################################
# Sensitivity and Specificity
##################################################

data_for_sensitivity_calc <- compare_results %>%
  filter(path_block_manual_check == "pathology blocks match")

true_positives <- nrow(data_for_sensitivity_calc %>%
                         filter(seqone_hrd_status == "POSITIVE" &
                                  hrd_status_check == "Seqone HRD status consistent with Myriad"))


true_negatives <- nrow(data_for_sensitivity_calc %>%
                         filter(seqone_hrd_status == "NEGATIVE" &
                                  hrd_status_check == "Seqone HRD status consistent with Myriad"))

false_positives <- nrow(data_for_sensitivity_calc %>%
                         filter(seqone_hrd_status == "POSITIVE" &
                                  hrd_status_check == "Seqone HRD status NOT consistent with Myriad"))

false_negatives <- nrow(data_for_sensitivity_calc %>%
                          filter(seqone_hrd_status == "NEGATIVE" &
                                   hrd_status_check == "Seqone HRD status NOT consistent with Myriad"))


data_table <- as.table(matrix(c(true_positives, false_positives, 
                                false_negatives, true_negatives), 
                              nrow = 2, byrow = TRUE))

data_results <- epiR::epi.tests(data_table, conf.level = 0.95)

##################################################
# QC metrics
##################################################

# Million reads
compare_results %>%
    filter(path_block_manual_check == "pathology blocks match") %>%
    ggplot(aes(x = reorder(shallow_sample_id, million_reads),
               y = million_reads)) +
    geom_point(size = 3, aes(colour = hrd_status_check)) +
  scale_colour_manual(values = c(negative_colour,
                                 positive_colour)) +
    labs(x = "", y = "million_reads") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "bottom",
        legend.title = element_blank())

# Read length
compare_results %>%
  filter(path_block_manual_check == "pathology blocks match") %>%
  ggplot(aes(x = reorder(shallow_sample_id, read_length),
             y = read_length)) +
  geom_point(size = 3, aes(colour = hrd_status_check)) +
  scale_colour_manual(values = c(negative_colour,
                                 positive_colour)) +
  labs(x = "", y = "read_length") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "bottom",
        legend.title = element_blank())

# Insert size

compare_results %>%
  filter(path_block_manual_check == "pathology blocks match") %>%
  ggplot(aes(x = reorder(shallow_sample_id, insert_size),
             y = insert_size)) +
  geom_point(size = 3, aes(colour = hrd_status_check)) +
  scale_colour_manual(values = c(negative_colour,
                                 positive_colour)) +
  labs(x = "", y = "insert_size") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "bottom",
        legend.title = element_blank())

# Percent Q30

compare_results %>%
  filter(path_block_manual_check == "pathology blocks match") %>%
  ggplot(aes(x = reorder(shallow_sample_id, percent_q30),
             y = percent_q30)) +
  geom_point(size = 3, aes(colour = hrd_status_check)) +
  scale_colour_manual(values = c(negative_colour,
                                 positive_colour)) +
  labs(x = "", y = "percent_q30") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "bottom",
        legend.title = element_blank())

# Percent aligned

compare_results %>%
  filter(path_block_manual_check == "pathology blocks match") %>%
  ggplot(aes(x = reorder(shallow_sample_id, percent_aligned),
             y = percent_aligned)) +
  geom_point(size = 3, aes(colour = hrd_status_check)) +
  scale_colour_manual(values = c(negative_colour,
                                 positive_colour)) +
  labs(x = "", y = "percent_aligned") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "bottom",
        legend.title = element_blank())

# Percent dups

compare_results %>%
  filter(path_block_manual_check == "pathology blocks match") %>%
  ggplot(aes(x = reorder(shallow_sample_id, percent_dups),
             y = percent_dups)) +
  geom_point(size = 3, aes(colour = hrd_status_check)) +
  scale_colour_manual(values = c(negative_colour,
                                 positive_colour)) +
  labs(x = "", y = "percent_dups") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "bottom",
        legend.title = element_blank())

# Coverage

compare_results %>%
  filter(path_block_manual_check == "pathology blocks match") %>%
  ggplot(aes(x = reorder(shallow_sample_id, coverage.y),
             y = coverage.y)) +
  geom_point(size = 3, aes(colour = hrd_status_check)) +
  scale_colour_manual(values = c(negative_colour,
                                 positive_colour)) +
  labs(x = "", y = "coverage.y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "bottom",
        legend.title = element_blank())

compare_results %>%
  filter(path_block_manual_check == "pathology blocks match") %>%
  ggplot(aes(x = insert_size,
             y = read_length)) +
  geom_point(size = 3, aes(colour = hrd_status_check)) +
  scale_colour_manual(values = c(negative_colour,
                                 positive_colour)) +
  labs(x = "insert_size", y = "read_length") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "bottom",
        legend.title = element_blank())

##################################################