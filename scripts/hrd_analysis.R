################################################################################
# HRD Analysis
# joseph.shaw3@nhs.net
################################################################################

##################################################
# Packages and Filepaths
##################################################

source("scripts/collate_seqone_pdfs.R")
source("scripts/collate_myriad_pdfs.R")

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

compare_results <- collated_myriad_info %>%
  left_join(seqone_mod, by = "nhs_number") %>%
  filter(!is.na(hrd_score))

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

comparison_plot <- ggplot(compare_results, aes(x = gi_score, y = hrd_score)) +
  geom_point(size = 2, alpha = 0.5) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(limits = c(0,100),
                     breaks = c(0, 25, 42, 50, 75, 100)) +
  ylim(0, 1) +
  labs(x = "Myriad Genome Instability Score",
       y = "SeqOne HRD Score",
       title = "Comparison of Myriad vs SeqOne HRD Testing") +
  geom_vline(xintercept = 42, linetype = "dashed") +
  geom_hline(yintercept = 0.5, linetype = "dashed")

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