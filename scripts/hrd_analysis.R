################################################################################
# HRD Analysis
# joseph.shaw3@nhs.net
################################################################################

##################################################
# Packages and Filepaths
##################################################

source("scripts/collate_seqone_pdfs.R")
source("scripts/collate_myriad_pdfs.R")

library("ggpubr")
library("readxl")

##################################################
# Additional data
##################################################

low_res_myriad_results <- read_excel(path = paste0(hrd_data_path,
                                                   "myriad_reports_low_res/myriad_low_res_report_data.xlsx")) %>%
  mutate(nhs_number = as.numeric(gsub(pattern = "(\\D)", "", nhs_no))) %>%
  select(-nhs_no)

collated_myriad_info_mod <- rbind(collated_myriad_info, low_res_myriad_results)

# Manual pathology ID check
path_block_check <- read_csv(paste0(hrd_project_path, "outputs/for_manual_check.csv")) %>%
  select(dlms_dna_number, check)

##################################################
# Compare SeqOne and Myriad Results
##################################################

dlms_204 <- read_csv(paste0(hrd_data_path, "DDBK_Samples_204.csv")) %>%
  janitor::clean_names() %>%
  dplyr::rename(dlms_dna_number = labno)

dlms_215 <- read_csv(paste0(hrd_data_path, "DDBK_Samples_215.csv")) %>%
  janitor::clean_names() %>%
  dplyr::rename(dlms_dna_number = labno) 

dlms_joined <- rbind(dlms_204, dlms_215)

seqone_mod <- collated_seqone_info %>%
  filter(!base::duplicated(dlms_dna_number)) %>%
  left_join(dlms_joined %>%
              select(dlms_dna_number, nhsno, firstname, surname, i_gene_r_no),
            by = "dlms_dna_number") %>%
  dplyr::rename(nhs_number = nhsno)

control_variants <- c("Biobank", "Seraseq")

control_names <- unique(grep(paste(control_variants,collapse="|"), seqone_mod$surname,
            value = TRUE, ignore.case = TRUE))

compare_results <- seqone_mod %>%
  left_join(collated_myriad_info_mod, by = "nhs_number") %>%
  left_join(path_block_check, by = "dlms_dna_number") %>%
  filter(!is.na(myriad_gi_score)) %>%
  mutate(check_outcome = case_when(
    
    myriad_hrd_status == "POSITIVE" & seqone_hrd_status == "POSITIVE" ~"consistent",
    myriad_hrd_status == "NEGATIVE" & seqone_hrd_status == "NEGATIVE" ~"consistent",
    myriad_hrd_status == "NEGATIVE" & seqone_hrd_status == "POSITIVE" ~"NOT consistent",
    myriad_hrd_status == "POSITIVE" & seqone_hrd_status == "NEGATIVE" ~"NOT consistent",
    TRUE ~"other"),
    
    identity = ifelse(surname %in% control_names, "Control", "Patient")) %>%
  filter(!identity == "Control" & dlms_dna_number != 21006723)

##################################################
# Comparison Plot
##################################################

# Run 1: 8 samples with Myriad results, 14 samples overall
# Run 2: 22 new samples with Myriad results, 31 samples overall
# 30 samples with Myriad results overall, 45 samples overall

ggplot(compare_results, aes(x = myriad_gi_score, y = seqone_hrd_score)) +
  geom_point(size = 2, aes(colour = check)) +
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
  facet_wrap(~check)

comparison_summary <- compare_results %>%
         group_by(check_outcome) %>%
         summarise(total = n())

ggplot(comparison_summary, aes(x = check_outcome, y = total)) +
  geom_col() +
  geom_text(aes(label = total), vjust = -0.5) +
  theme_bw() +
  labs(y = "Total Samples", x = "", title = "Consistency of Myriad and Seqone Testing")

##################################################
# Repeat Testing
##################################################

collated_seqone_info %>%
  filter(base::duplicated(dlms_dna_number, fromLast = TRUE) |
           base::duplicated(dlms_dna_number, fromLast = FALSE)) %>%
  ggplot(aes(x = filename, y = hrd_score)) +
  geom_point() +
  facet_wrap(~dlms_dna_number) +
  theme(axis.text.x = element_blank())

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

ggplot(seq_one_results, aes(x = reorder(filename, hrd_score),
                            y = hrd_score)) +
  geom_point(size = 2, aes(colour = result_model)) +
  scale_colour_manual(values = c(negative_colour,
                                 positive_colour)) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "", y = "SeqOne HRD Score",
       title = "Explaining the SeqOne HRD Model")

ggplot(seq_one_results, aes(x = lga,
                            y = lpc,
                            colour = hrd_status)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(x = "Large Genomic Alterations", 
       y = "Loss of Parental Copy") +
  geom_vline(xintercept = 18, linetype = "dashed") +
  geom_vline(xintercept = 14, linetype = "dashed") +
  geom_hline(yintercept = 10, linetype = "dashed")

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