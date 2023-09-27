################################################################################
# HRD Analysis
# joseph.shaw3@nhs.net
################################################################################

rm(list=ls())

##################################################
# Packages and Filepaths
##################################################

library("ggpubr")
library("readxl")

##################################################
# Source scripts
##################################################

source("scripts/collate_seqone_pdfs.R")
source("scripts/collate_myriad_pdfs.R")
source("scripts/dlms_connection.R")

downsampled_samples <- grep(pattern = "downsampl", x = collated_seqone_info$sample_id, value = TRUE)

##################################################
# Collate and edit Myriad Information
##################################################

# pdf_text doesn't work on these reports as resolution is too low.
# I had to collate data manually in an Excel
low_res_myriad_results <- read_excel(path = paste0(hrd_data_path,
                                                   "myriad_reports_low_res/myriad_low_res_report_data.xlsx")) %>%
  mutate(nhs_number = as.numeric(gsub(pattern = "(\\D)", "", nhs_number))) 

collated_myriad_info_mod <- rbind(collated_myriad_info, low_res_myriad_results)

##################################################
# Pathology Block ID Check
##################################################

# Extract data from DLMS via ODBC
seqone_dlms_info <- get_sample_data(collated_seqone_info$dlms_dna_number) %>%
  dplyr::rename(dlms_dna_number = labno,
                nhs_number = nhsno) 

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

write.csv(export_for_check, paste0(hrd_data_path, "manual_path_block_check.csv"),
          row.names = FALSE)


path_block_check <- read_csv(paste0(hrd_data_path, "manual_path_block_check.csv")) %>%
  select(dlms_dna_number, path_block_manual_check)

##################################################
# Compare SeqOne and Myriad Results
##################################################

seqone_mod <- collated_seqone_info %>%
  filter(!date %in% c("September 1, 2023", "August 31, 2023",
                      "August 24, 2023", "August 25, 2023")) %>%
  left_join(seqone_dlms_info %>%
              select(dlms_dna_number, nhs_number, firstname, surname, i_gene_r_no,
                     pathno),
            by = "dlms_dna_number") %>%
  mutate(downsampled = ifelse(sample_id %in% downsampled_samples, "Yes", "No"))

control_variants <- c("Biobank", "Seraseq")

control_names <- unique(grep(paste(control_variants,collapse="|"), seqone_mod$surname,
            value = TRUE, ignore.case = TRUE))

compare_results <- seqone_mod %>%
  filter(!base::duplicated(dlms_dna_number)) %>%
  left_join(collated_myriad_info_mod, by = "nhs_number") %>%
  left_join(path_block_check, by = "dlms_dna_number") %>%
  #filter(!is.na(myriad_gi_score)) %>%
  mutate(hrd_status_check = case_when(
    
    myriad_hrd_status == "POSITIVE" & seqone_hrd_status == "POSITIVE" ~"consistent",
    myriad_hrd_status == "NEGATIVE" & seqone_hrd_status == "NEGATIVE" ~"consistent",
    myriad_hrd_status == "NEGATIVE" & seqone_hrd_status == "POSITIVE" ~"NOT consistent",
    myriad_hrd_status == "POSITIVE" & seqone_hrd_status == "NEGATIVE" ~"NOT consistent",
    TRUE ~"other"),
    
    identity = ifelse(surname %in% control_names, "Control", "Patient"))

  # Rearrange columns for easier viewing
  select(dlms_dna_number, sample_id, i_gene_r_no, myriad_r_number,
         nhs_number, 
         firstname, surname, myriad_patient_name,
         seqone_hrd_score, myriad_gi_score, seqone_hrd_status,
         myriad_hrd_status, hrd_status_check,
         pathno, myriad_pathology_block, pathology_block_check)

# filter(!identity == "Control" & dlms_dna_number != 21006723)

##################################################
# Comparison Plot
##################################################

# Run 1: 8 samples with Myriad results, 14 samples overall
# Run 2: 22 new samples with Myriad results, 31 samples overall
# 30 samples with Myriad results overall, 45 samples overall

ggplot(compare_results %>%
         filter(!is.na(path_block_manual_check)), aes(x = myriad_gi_score, y = seqone_hrd_score)) +
  geom_point(size = 2, aes(shape = hrd_status_check)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(limits = c(0,100),
                     breaks = c(0, 25, 42, 50, 75, 100)) +
  ylim(0, 1) +
  labs(x = "Myriad Genome Instability Score",
       y = "SeqOne HRD Score",
       title = "Comparison of Myriad vs SeqOne HRD Testing",
       #subtitle = paste0(length(unique(compare_results$dlms_dna_number)), " samples")
       ) +
  #geom_vline(xintercept = 42, linetype = "dashed") +
  #geom_hline(yintercept = 0.5, linetype = "dashed") + 
  ggpubr::stat_cor(method = "pearson", label.x = 50, label.y = 0.25) +
  facet_wrap(~path_block_manual_check)

comparison_summary <- compare_results %>%
  filter(!is.na(myriad_gi_score)) %>%
   group_by(hrd_status_check) %>%
   summarise(total = n())

ggplot(comparison_summary, aes(x = hrd_status_check, y = total)) +
  geom_col() +
  geom_text(aes(label = total), vjust = -0.5) +
  theme_bw() +
  labs(y = "Total Samples", x = "", title = "Consistency of Myriad and Seqone Testing")

inconsistent_summary <- compare_results %>%
  filter(hrd_status_check =="NOT consistent") %>%
  select(dlms_dna_number, myriad_patient_name, seqone_hrd_score, myriad_gi_score,
         seqone_hrd_status, myriad_hrd_status, pathno, myriad_pathology_block, 
         pathology_block_check) %>%
  arrange(pathology_block_check)

write.csv(inconsistent_summary,
          paste0(hrd_project_path, "outputs/inconsistent_summary.csv"),
          row.names = FALSE)


##################################################
# Repeat Testing
##################################################

seqone_mod %>%
  filter(base::duplicated(dlms_dna_number, fromLast = TRUE) |
           base::duplicated(dlms_dna_number, fromLast = FALSE)) %>%
  ggplot(aes(x = worksheet, y = seqone_hrd_score)) +
  geom_jitter(size = 3, alpha = 0.5, aes(colour = downsampled)) +
  facet_wrap(~dlms_dna_number) +
  theme_bw() +
  #theme(axis.text.x = element_blank()) +
  labs(title = "SeqOne results for repeated samples",
       x = "Worksheet",
       y = "SeqOne HRD score") +
  geom_hline(yintercept = 0.50, linetype = "dashed")

unusual_samples <- c(21013520, 23032088, 20127786)

unusual_repeat_sample_info <- seqone_mod %>%
  filter(dlms_dna_number %in% unusual_samples) %>%
  select(-c(sample_id, filename, date, pathno)) %>%
  arrange(dlms_dna_number)

write.csv(unusual_repeat_sample_info,
          paste0(hrd_project_path, "outputs/unusual_repeat_sample_info.csv"),
          row.names = FALSE)

##################################################
# Comparison with Simplified Model
##################################################

seq_one_results <- seqone_mod %>%
  mutate(result_model = case_when(
    lga >= 18 ~"Positive",
    lga <= 14 ~"Negative",
    lga < 18 & lga > 14 & lpc > 10 ~"Positive",
    lga < 18 & lga > 14 & lpc <= 10 ~"Negative")) %>%
  filter(!base::duplicated(dlms_dna_number))

##################################################
# Plots
##################################################

positive_colour <- "#FF3333"
negative_colour <- "#3300FF"

ggplot(seq_one_results, aes(x = reorder(filename, seqone_hrd_score),
                            y = seqone_hrd_score)) +
  geom_point(size = 2, aes(colour = result_model)) +
  scale_colour_manual(values = c(negative_colour,
                                 positive_colour)) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  labs(x = "", y = "SeqOne HRD Score",
       title = "Explaining the SeqOne HRD Model")

ggplot(seq_one_results, aes(x = lga,
                            y = lpc,
                            colour = seqone_hrd_status)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(x = "Large Genomic Alterations", 
       y = "Loss of Parental Copy") +
  geom_vline(xintercept = 18, linetype = "dashed") +
  geom_vline(xintercept = 14, linetype = "dashed")
  # geom_hline(yintercept = 10, linetype = "dashed")

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
# Run 1
##################################################

WS133557 <- seqone_mod %>%
  filter(worksheet == "WS133557")


total_samples <- length(unique(WS133557$dlms_dna_number))

# Total places on the run = 20
# 3 samples repeated 3 times = 9
# (14-3) + 9 = 20

total_analyses <- nrow(WS133557)

down_sampling <- collated_seqone_info %>%
  filter(worksheet == "WS133557" & sample_id %in% downsampled_samples)

collated_seqone_info %>%
  filter(worksheet == "WS133557") %>%
  group_by(date) %>%
  summarise(total = n()) %>%
  arrange(total)

check <- collated_seqone_info %>%
  filter(date %in% c("September 6, 2023"))

##################################################
# Seraseq Controls
##################################################

seraseq_gi_scores <- read_excel(paste0(hrd_data_path, "seraseq_gi_scores.xlsx"))

seraseq_control_data <- seqone_mod %>%
  left_join(seraseq_gi_scores, by = "dlms_dna_number") %>%
  filter(!is.na(myriad_gi_score)) %>%
  mutate(firstname_factor = factor(firstname, levels = c(
    "FFPE HRD Negative",
    "Low-Positive FFPE HRD",
    "High-Positive FFPE HRD")))

ggplot(seraseq_control_data, aes(x = myriad_gi_score, y = seqone_hrd_score)) +
  geom_point(size = 3, aes(shape = seqone_hrd_status)) +
  ylim(0, 1) +
  xlim(0, 100) +
  facet_wrap(~firstname_factor) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "Myriad GI Score", y = "SeqOne HRD Score")

ggplot(seraseq_control_data, aes(x = lga, y = percent_mapping)) +
  geom_point(size = 3, aes(shape = seqone_hrd_status)) +
  facet_wrap(~firstname_factor)


##################################################