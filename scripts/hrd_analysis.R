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
source("functions/hrd_functions.R")

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
         myriad_pathology_block_pg1 = "",
         myriad_pathology_block_pg2 = "",
         myriad_brca_status = "",
         myriad_filename = "") %>%
  select(myriad_r_number, myriad_patient_name, myriad_dob, nhs_number,
         myriad_pathology_block_pg1, myriad_pathology_block_pg2,
         myriad_gi_score, myriad_hrd_status, 
         myriad_brca_status, myriad_filename)

collated_myriad_info_mod <- rbind(collated_myriad_info, low_res_myriad_results,
                                  seraseq_gi_scores)

##################################################
# Pathology Block ID Check
##################################################

comment_regex_single <- ".+(\\W{1}\\d{2}%).+"

comment_regex_range <- ".+(\\d{2}\\W{1}\\d{2}%).+"

# Extract data from DLMS via ODBC
seqone_dlms_info <- get_sample_data(collated_seqone_info$dlms_dna_number) %>%
  dplyr::rename(dlms_dna_number = labno,
                nhs_number = nhsno) %>%
  mutate(ncc_single = sub(x = comments,
                          pattern = comment_regex_single,
                          replacement = "\\1"),
         ncc_range = sub(x = comments,
                         pattern = comment_regex_range,
                         replacement = "\\1"),
         ncc = ifelse(str_length(ncc_single) == 4 & ncc_range > 6,
                      ncc_single, ifelse(str_length(ncc_range) == 6, 
                                         ncc_range, NA)))

# Enter a fake NHS number for the Seraseq and Biobank controls, to allow
# Myriad scores to be joined later, and to keep Biobank controls within the 
# dataset

seqone_dlms_info[seqone_dlms_info$dlms_dna_number == 23032086, "nhs_number"] <- 1
seqone_dlms_info[seqone_dlms_info$dlms_dna_number == 23032088, "nhs_number"] <- 2
seqone_dlms_info[seqone_dlms_info$dlms_dna_number == 23031639, "nhs_number"] <- 3
seqone_dlms_info[seqone_dlms_info$dlms_dna_number == 23033285, "nhs_number"] <- 4
seqone_dlms_info[seqone_dlms_info$dlms_dna_number == 23033279, "nhs_number"] <- 5
seqone_dlms_info[seqone_dlms_info$dlms_dna_number == 23033288, "nhs_number"] <- 6

export_for_check <- collated_seqone_info %>%
  filter(!base::duplicated(dlms_dna_number)) %>%
  select(dlms_dna_number) %>%
  left_join(seqone_dlms_info %>%
              select(dlms_dna_number, firstname, surname, pathno, nhs_number), 
            by = "dlms_dna_number") %>%
  left_join(collated_myriad_info_mod %>%
              select(nhs_number, myriad_pathology_block_pg1, 
                     myriad_pathology_block_pg2), by = "nhs_number") %>%
  select(dlms_dna_number, firstname, surname, nhs_number, 
         pathno, myriad_pathology_block_pg1, myriad_pathology_block_pg2) %>%
  arrange(nhs_number) %>%
  mutate(path_block_automated_check = ifelse(pathno == myriad_pathology_block_pg2,
                                             "match", "NOT match"),
         path_block_manual_check = "")

# Writing of pathology block IDs is inconsistent, so a manual check is
# required

export_timestamp(hrd_data_path, export_for_check)

path_block_check <- read_excel(paste0(hrd_data_path, "manual_path_block_check_edit.xlsx")) %>%
  select(dlms_dna_number, path_block_manual_check)

##################################################
# Initial DNA Concentrations
##################################################

dna_concentrations <- read_excel(path = paste0(hrd_data_path, "HS2_sample_prep_export.xlsx"),
                                 col_types = c("date", "numeric","text", "text",
                                               "date", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric",
                                               "date", "text", "numeric",
                                               "numeric", "text", "text")) %>%
  janitor::clean_names() %>%
  dplyr::rename(dlms_dna_number = sample_id,
                forename = name,
                surname = x4,
                nanodrop_ng_ul = nanodrop,
                nd_260_280 = x7,
                nd_260_230 = x8,
                stock_volume_ul = qubit,
                qubit_dna_ul = x10,
                dilution_date = dilution,
                dilution_worksheet = x12,
                dilution_dna_volume = x13,
                dilution_water_volume = x14) %>%
  filter(!is.na(dlms_dna_number))

neat_variants <- unique(grep(pattern = "neat", 
                             dna_concentrations$comments, value = TRUE))

dna_concentrations_mod <- dna_concentrations %>%

  mutate(dilution_concentration = ifelse(
    
    comments %in% neat_variants,
    
    qubit_dna_ul,
    
    round((dilution_dna_volume * qubit_dna_ul) / (dilution_water_volume + dilution_dna_volume), 2)),
    
    # 15ul of diluted or neat (undiluted) DNA used in fragmentation reaction
    input_ng = dilution_concentration * 15,
    # 1 haploid genome = 3.3 picograms
    input_genomes = (input_ng*1000)/3.3) %>%
  # Some samples used on previous runs. Filtering by comments isolates
  # the dilutions for the HRD runs
  filter(!is.na(comments)) %>%
  # Checked: repeat rows of the same DNA number have the same DNA concentration
  # and dilution volumes
  filter(!base::duplicated(dlms_dna_number))

##################################################
# qPCR Library QC
##################################################

hs2_library_prep <- read_excel(path = "S:/central shared/Genetics/Repository/Technical Teams/NGS/SureSelect XT HS2/SSXT HS2 Library Prep 2023.xlsx",
                               sheet = "Sheet1",
                               col_types = c("text", "guess", "text", "text",
                                             "text", "numeric", "numeric",
                                             "numeric", "numeric", "numeric", 
                                             "numeric", "numeric", "numeric",
                                             "numeric", "numeric", "numeric",
                                             "numeric", "numeric", "text",
                                             "guess", "guess", "guess", "guess",
                                             "guess", "guess", "guess", "guess")) %>%
  janitor::clean_names() %>%
  filter(!plate_position %in% c("Plate Position", NA))

kapa_data_WS133557 <- extract_kapa_data("133557", 20)

kapa_data_WS134687 <- extract_kapa_data("134687", 31)

kapa_data_WS134928 <- extract_kapa_data("134928", 7)

kapa_data_WS135001 <- extract_kapa_data("135001", 31)

kapa_data_collated <- rbind(kapa_data_WS133557, kapa_data_WS134687,
                            kapa_data_WS134928, kapa_data_WS135001)

##################################################
# Compare SeqOne and Myriad Results
##################################################

# Add QC data for SeqOne samples

seqone_qc_data <- read_excel(paste0(hrd_data_path, "seqone_qc_metrics/seqone_qc_metrics_2023_10_03.xlsx")) %>%
  janitor::clean_names() %>%
  dplyr::rename(shallow_sample_id = sample,
                read_length = read_len,
                insert_size = ins_size,
                million_reads = m_reads)

seqone_mod <- collated_seqone_info %>%
  filter(!date %in% c("September 1, 2023", "August 31, 2023",
                      "August 24, 2023", "August 25, 2023")) %>%
  left_join(seqone_dlms_info %>%
              select(dlms_dna_number, nhs_number, firstname, surname,
                     i_gene_r_no, pathno, ncc),
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
  filter(downsampled == "No") %>%
  left_join(dna_concentrations_mod %>%
              select(dlms_dna_number, qubit_dna_ul, input_ng,
                     input_genomes), by = "dlms_dna_number") %>%
  left_join(kapa_data_collated %>%
              select(shallow_sample_id, q_pcr_n_m, ts_ng_ul, total_yield, q_pcr_n_m),
            by = "shallow_sample_id")

# Classify Seraseq controls as "pathology blocks match"
compare_results[compare_results$dlms_dna_number == 23032086, "path_block_manual_check"] <- "pathology blocks match"
compare_results[compare_results$dlms_dna_number == 23032088, "path_block_manual_check"] <- "pathology blocks match"
compare_results[compare_results$dlms_dna_number == 23031639, "path_block_manual_check"] <- "pathology blocks match"

##################################################
# Repeat Testing Plots - Inter Run Variation
##################################################

repeat_results <- compare_results %>%
  filter(base::duplicated(dlms_dna_number, fromLast = TRUE) |
           base::duplicated(dlms_dna_number, fromLast = FALSE)) %>%
  filter(downsampled == "No") %>%
  mutate(input_category = case_when(
    input_ng >= 48 ~"50ng input",
    input_ng < 48 ~"lower than 50ng input"))


repeat_facet_plot <- ggplot(repeat_results, aes(x = worksheet, 
                                                y = seqone_hrd_score)) +
                      geom_point(alpha = 0.5, 
                                 aes(shape = seqone_hrd_status,
                                     size = coverage.x)) +
                      facet_wrap(~dlms_dna_number) +
                      theme_bw() +
                      theme(axis.text.x = element_text(angle = 90),
                            legend.position = "bottom") +
                      labs(title = "SeqOne results for repeated samples",
                           x = "",
                           y = "SeqOne HRD score",
                           caption = "Plot name: repeat_facet_plot") +
                      geom_hline(yintercept = 0.50, linetype = "dashed")

save_hrd_plot(repeat_facet_plot, input_width = 15, input_height = 15)

sample_20127786_plot <- make_individual_plot(20127786)

sample_21011999_plot <- make_individual_plot(21011999)

sample_23032088_plot <- make_individual_plot(23032088)

sample_21013520_plot <- make_individual_plot(21013520)

##################################################
# Intra Run Variation
##################################################

intra_run_plot <- compare_results %>%
  filter(worksheet == "WS133557") %>%
  filter(base::duplicated(dlms_dna_number, fromLast = TRUE) |
           base::duplicated(dlms_dna_number, fromLast = FALSE)) %>%
  ggplot(aes(x = sample_id,
             y = seqone_hrd_score)) +
           geom_point(size = 3) +
  labs(x = "", y = "SeqOne HRD Score",
       title = "Repeated Sample Results",
       subtitle = "Worksheet WS133557") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))

save_hrd_plot(intra_run_plot, input_width = 10, input_height = 10)

##################################################
# Repeat Testing
##################################################

summary_21013520 <- get_sample_summary_info(21013520)

export_timestamp(hrd_output_path, summary_21013520)

summary_23032088 <- get_sample_summary_info(23032088)

export_timestamp(hrd_output_path, summary_23032088)

summary_20127786 <- get_sample_summary_info(20127786)

export_timestamp(hrd_output_path, summary_20127786)

summary_21011999 <- get_sample_summary_info(21011999)

export_timestamp(hrd_output_path, summary_21011999)

##################################################
# Compare Results of Mid Output Run - WS134928
##################################################

mid_output_run <- seqone_mod %>%
  filter(worksheet == "WS134928")

WS134928_repeated_samples <- mid_output_run$dlms_dna_number

WS134928_WS134687_data <- compare_results %>%
  filter(dlms_dna_number %in% WS134928_repeated_samples &
           worksheet %in% c("WS134687", "WS134928"))

WS134928_WS134687_plot <- ggplot(WS134928_WS134687_data, 
                                 aes(x = coverage.x, y = seqone_hrd_score)) +
  geom_point(size = 3, aes(shape = seqone_hrd_status,
                           colour = worksheet)) +
  facet_wrap(~dlms_dna_number) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "Coverage", y = "SeqOne HRD Score",
       title = "Sample libraries repeated on mid-output run (WS134928)",
       caption = "Plot name: WS134928_WS134687_plot") +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  xlim(0.8, 1.8)

##################################################
# HRD Status Comparison with Myriad
##################################################

comparison_summary <- compare_results %>%
  filter(!is.na(myriad_gi_score) & myriad_gi_score != "NA") %>%
  filter(path_block_manual_check != "NA") %>%
  group_by(hrd_status_check, path_block_manual_check) %>%
  summarise(total = n())

myriad_comparison_plot <- ggplot(comparison_summary, aes(x = hrd_status_check, 
                                                         y = total)) +
  geom_col() +
  geom_text(aes(label = total), vjust = -0.5) +
  theme_bw() +
  labs(y = "Total DNA Inputs", x = "", 
       title = "Consistency of Myriad and Seqone Testing",
       subtitle = paste0("Data for ", sum(comparison_summary$total), " DNA inputs"),
       caption = "Plot name: myriad_comparison_plot") +
  facet_wrap(~path_block_manual_check) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##################################################
# Impact of Pathology Blocks
##################################################

results_for_path_block_plot <- compare_results %>%
  filter(path_block_manual_check != "NA" & !is.na(myriad_gi_score))

path_block_plot <- ggplot(results_for_path_block_plot, 
                          aes(x = myriad_gi_score, y = seqone_hrd_score)) +
  geom_point(size = 3, alpha = 0.6, aes(shape = hrd_status_check)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(limits = c(0,100),
                     breaks = c(0, 25, 42, 50, 75, 100)) +
  ylim(0, 1) +
  labs(x = "Myriad Genome Instability Score",
       y = "SeqOne HRD Score",
       title = "Comparison of Myriad vs SeqOne HRD Testing",
       subtitle = paste0("Data for ", nrow(results_for_path_block_plot), " DNA inputs")) +
  ggpubr::stat_cor(method = "pearson", label.x = 50, label.y = 0.25) +
  #facet_wrap(~path_block_manual_check) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom",
        legend.title = element_blank()) +
  guides(shape=guide_legend(ncol=1)) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 42, linetype = "dashed")

high_quality_results <- compare_results %>%
  filter(path_block_manual_check != "NA" & !is.na(myriad_gi_score)) %>%
  filter(input_ng >= 49 & coverage.x >= 1 &
           path_block_manual_check == "pathology blocks match")

# Results with quality control added
ggplot(high_quality_results, 
       aes(x = myriad_gi_score, y = seqone_hrd_score)) +
  geom_point(size = 3, alpha = 0.6, aes(shape = hrd_status_check)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(limits = c(0,100),
                     breaks = c(0, 25, 42, 50, 75, 100)) +
  ylim(0, 1) +
  labs(x = "Myriad Genome Instability Score",
       y = "SeqOne HRD Score",
       title = "Comparison of Myriad vs SeqOne HRD Testing",
       subtitle = paste0("Data for ", nrow(high_quality_results), " DNA inputs: path blocks match, >=49ng input, >=1x coverage")) +
  ggpubr::stat_cor(method = "pearson", label.x = 50, label.y = 0.25) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom",
        legend.title = element_blank()) +
  guides(shape=guide_legend(ncol=1)) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 42, linetype = "dashed")

# Coverage vs DNA input
ggplot(compare_results %>%
         filter(path_block_manual_check == "pathology blocks match"), 
       aes(x = input_ng, y = coverage.x)) +
  geom_point(size = 3, alpha = 0.6, aes(shape = hrd_status_check,
                                        colour = hrd_status_check)) +
  scale_colour_manual(values = c("#CCCCCC", "#FF0000")) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom",
        legend.title = element_blank()) +
  guides(shape=guide_legend(ncol=1)) +
  ylim(0, 2)

# Total yield vs DNA input
ggplot(compare_results %>%
         filter(path_block_manual_check == "pathology blocks match"), 
       aes(x = input_ng, y = total_yield)) +
  geom_point(size = 3, alpha = 0.6, aes(shape = hrd_status_check,
                                        colour = hrd_status_check)) +
  scale_colour_manual(values = c("#CCCCCC", "#FF0000")) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom",
        legend.title = element_blank())

# DNA samples with lower than 50ng input
ggplot(compare_results %>%
         filter(input_ng < 49 & !is.na(myriad_gi_score)), 
       aes(x = myriad_gi_score, y = seqone_hrd_score)) +
  geom_point(size = 3, alpha = 0.6, aes(shape = hrd_status_check,
                                        colour = hrd_status_check)) +
  scale_colour_manual(values = c("#CCCCCC", "#FF0000")) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom",
        legend.title = element_blank()) +
  xlim(0, 100) +
  ylim(0, 1) +
  ggpubr::stat_cor(method = "pearson", label.x = 50, label.y = 0.25) +
  labs(title = "DNA samples with lower than 50ng input")

##################################################
# Individual discrepant samples
##################################################

# Sample 23034142
red_circle_23034142 <- circle_individual_point(23034142)

summary_23034142 <- get_sample_summary_info(23034142)

export_timestamp(hrd_output_path, summary_23034142)

# Sample 23016518
red_circle_23016518 <- circle_individual_point(23016518)

summary_23016518 <- get_sample_summary_info(23016518)

export_timestamp(hrd_output_path, summary_23016518)

# Sample 23016516
red_circle_23016516 <- circle_individual_point(23016516)

summary_23016516 <- get_sample_summary_info(23016516)

export_timestamp(hrd_output_path, summary_23016516)

# Sample 23016526
red_circle_23016526 <- circle_individual_point(23016526)

summary_23016526 <- get_sample_summary_info(23016526)

export_timestamp(hrd_output_path, summary_23016526)

coverage_plot <- compare_results %>%
  filter(path_block_manual_check != "NA") %>%
  ggplot(aes(x = reorder(shallow_sample_id, coverage.x), y = coverage.x)) +
  geom_point(size = 3, aes(shape = hrd_status_check,
                           colour = hrd_status_check)) +
  scale_colour_manual(values = c("#CCCCCC", "#FF0000")) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank()) +
  labs(x = "Sample", y = "Coverage",
       caption = "Plot name: coverage_plot") +
  facet_wrap(~path_block_manual_check)

##################################################
# Export tables
##################################################

inconsistent_summary <- compare_results %>%
  filter(hrd_status_check == "Seqone HRD status NOT consistent with Myriad") %>%
  select(sample_id, worksheet, seqone_hrd_score, myriad_gi_score,
         seqone_hrd_status, myriad_hrd_status,
         path_block_manual_check, lga, lpc, ccne1, rad51b, coverage.x, percent_mapping,
         myriad_r_number, input_ng) %>%
  arrange(path_block_manual_check, sample_id)

export_timestamp(hrd_output_path, inconsistent_summary)

# Entire table of results
export_timestamp(hrd_output_path, compare_results)

# Remove patient identifiable information

results_patient_ids_removed <- compare_results %>%
  select(-c(nhs_number, firstname, surname, myriad_patient_name, 
            myriad_dob, myriad_filename))

export_timestamp(hrd_output_path, results_patient_ids_removed)

##################################################
# Impact of read length
##################################################

compare_results %>%
  filter(downsampled == "No" & path_block_manual_check != "NA") %>%
  ggplot(aes(x = reorder(shallow_sample_id, read_length), y = read_length)) +
  geom_point(size = 3, alpha = 0.5, aes(shape = hrd_status_check)) +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  labs(title = "SeqOne Read Length",
       x = "",
       y = "Read length") +
  facet_wrap(~path_block_manual_check)

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
                            y = lpc)) +
  #scale_fill_manual(values = c(negative_colour,
                               #positive_colour)) +
  geom_point(size = 3, alpha = 0.6,
             aes(colour = seqone_hrd_status)) +
  scale_shape_manual(values = c(24, 21)) +
  theme_bw() +
  labs(x = "Large Genomic Alterations", 
       y = "Loss of Parental Copy",
       title = "Comparison of Seqone model approximation to pipeline output",
       subtitle = paste0("Data for ", nrow(simplified_model), " BRCA-negative DNA inputs")) +
  geom_hline(yintercept = 10, linetype= "dashed") +
  geom_vline(xintercept = 20, linetype = "dashed")

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
# Biobank Controls
##################################################

compare_results %>%
  filter(surname %in% biobank_names) %>%
  ggplot(aes(x = surname, y = seqone_hrd_score)) +
  geom_point(size = 3) +
  ylim(0, 1) +
  labs(x = "", y = "SeqOne HRD score") +
  theme_bw()

##################################################
# Myriad results
##################################################

collated_myriad_info_mod %>%
  ggplot(aes(x = reorder(nhs_number, myriad_gi_score),
             y = myriad_gi_score)) +
  geom_point(size = 3, aes(shape = myriad_hrd_status)) +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  geom_hline(yintercept = 42, linetype = "dashed")

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
ggplot(compare_results, aes(x = reorder(shallow_sample_id, million_reads),
               y = million_reads)) +
    geom_point(size = 3, aes(colour = hrd_status_check)) +
    labs(x = "", y = "million_reads") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "bottom",
        legend.title = element_blank())

# Read length
ggplot(compare_results, aes(x = reorder(shallow_sample_id, read_length),
             y = read_length)) +
  geom_point(size = 3, aes(colour = hrd_status_check)) +
  #scale_colour_manual(values = c(negative_colour,
                                 #positive_colour)) +
  labs(x = "", y = "read_length") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "bottom",
        legend.title = element_blank()) +
  facet_wrap(~worksheet)

# Insert size
ggplot(compare_results, aes(x = reorder(shallow_sample_id, insert_size),
             y = insert_size)) +
  geom_point(size = 3, aes(colour = hrd_status_check)) +
  labs(x = "", y = "insert_size") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "bottom",
        legend.title = element_blank())

# Percent Q30
ggplot(compare_results, aes(x = reorder(shallow_sample_id, percent_q30),
             y = percent_q30)) +
  geom_point(size = 3, aes(colour = hrd_status_check)) +
  labs(x = "", y = "percent_q30") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "bottom",
        legend.title = element_blank()) +
  ylim(0, 100)

# Percent aligned
ggplot(compare_results, aes(x = reorder(shallow_sample_id, percent_aligned),
             y = percent_aligned)) +
  geom_point(size = 3, aes(colour = hrd_status_check)) +
  labs(x = "", y = "percent_aligned") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "bottom",
        legend.title = element_blank())

# Percent dups
ggplot(compare_results, aes(x = reorder(shallow_sample_id, percent_dups),
             y = percent_dups)) +
  geom_point(size = 3, aes(colour = hrd_status_check)) +
  labs(x = "", y = "percent_dups") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "bottom",
        legend.title = element_blank())

# Coverage
ggplot(compare_results, aes(x = reorder(shallow_sample_id, coverage.y),
             y = coverage.y)) +
  geom_point(size = 3, aes(colour = hrd_status_check)) +
  labs(x = "", y = "coverage.y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "bottom",
        legend.title = element_blank())

ggplot(compare_results, aes(x = coverage.x,
             y = read_length)) +
  geom_point(size = 3, aes(colour = hrd_status_check)) +
  labs(x = "coverage", y = "read_length") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "bottom",
        legend.title = element_blank())

##################################################
# Sample Extraction Information
##################################################

# Get tissues
seqone_tissue_types <- unique(seqone_dlms_info$tissue)

tissue_query <- sqlQuery(channel = moldb_connection,
                         query = paste0("SELECT * FROM MolecularDB.dbo.TissueTypes WHERE TissueTypeId IN (",
                                        paste(seqone_tissue_types, collapse = ", "),
                                        ")")) %>%
  janitor::clean_names()

# Get extraction batch IDs
extraction_batch_id_table <- sqlQuery(channel = moldb_connection,
                                  query = paste0("SELECT * FROM MolecularDB.dbo.MOL_Extractions WHERE LABNO IN (",
                                                 paste(seqone_dlms_info$dlms_dna_number, collapse = ", "),
                                                 ")")) %>%
  arrange(LabNo) %>%
  janitor::clean_names() %>%
  dplyr::rename(dlms_dna_number = lab_no,
                extraction_batch_id = extraction_batch_fk)


extraction_batches <- unique(extraction_batch_id_table$extraction_batch_id)

# Get extraction batch dates
extraction_batch_date_table <- sqlQuery(channel = moldb_connection,
                                   query = paste0("SELECT * FROM MolecularDB.dbo.MOL_ExtractionBatches WHERE ExtractionBatchId IN (",
                                                  paste(extraction_batches, collapse = ", "),
                                                  ")")) %>%
  janitor::clean_names()

extraction_batch_info <- extraction_batch_id_table %>%
  left_join(extraction_batch_date_table, by = "extraction_batch_id") %>%
  # Get only extraction batch IDs used for Cobas extractions
  filter(extraction_method_fk == 25)

# 1 DNA number (23033279) is on 2 extraction batches

seqone_dlms_extractions <- seqone_dlms_info %>%
  left_join(tissue_query %>%
              dplyr::rename(tissue = tissue_type_id), by = "tissue") %>%
  left_join(extraction_batch_info, by = "dlms_dna_number")

# Plot of tissue types
ggplot(seqone_dlms_extractions, aes(x = tissue_type, y = )) +
  geom_bar()

ggplot(seqone_dlms_extractions, aes(x = run_date, y = concentration)) +
  geom_point()

##################################################
# Impact of DNA Concentration
##################################################

ggplot(compare_results %>%
         filter(path_block_manual_check == "pathology blocks match"), 
       aes(x = input_ng, y = coverage.x)) +
  geom_point(size = 3, aes(shape = hrd_status_check,
                 colour = hrd_status_check)) +
  scale_colour_manual(values = c("#CCCCCC", "#FF0000")) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  ylim(0,2) +
  xlim(0, 55) +
  facet_wrap(~path_block_manual_check) 
  #geom_point(data=compare_results[compare_results$dlms_dna_number == 23034142,], 
             #aes(input_ng, coverage.x),
             #=21, fill=NA, size=5, colour="red", stroke=2)

ggplot(compare_results %>%
         filter(path_block_manual_check == "pathology blocks match"), 
       aes(x = input_genomes, y = lga)) +
  geom_point(size = 3, aes(shape = hrd_status_check,
                           colour = hrd_status_check)) +
  scale_colour_manual(values = c("#CCCCCC", "#FF0000")) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  facet_wrap(~path_block_manual_check)

ggplot(compare_results %>%
         filter(path_block_manual_check == "pathology blocks match"), 
       aes(x = reorder(dlms_dna_number, qubit_dna_ul),
                               y = qubit_dna_ul)) +
  geom_point(size = 3, aes(shape = hrd_status_check,
                           colour = hrd_status_check)) +
  scale_colour_manual(values = c("#CCCCCC", "#FF0000")) +
  geom_hline(yintercept = 3.3, linetype = "dashed") +
  scale_y_continuous(trans='log10') +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  labs(x = "") 

##################################################
# Tumour BRCA Referral DNA Concentrations
##################################################

brca_query <- "SELECT * FROM MolecularDB.dbo.Samples WHERE DISEASE IN (204)"

tbrca_data <- sqlQuery(channel = moldb_connection,
                        query = brca_query) %>%
  janitor::clean_names()

dna_qc_threshold <- round(100/15, 0)

tbrca_data_mod <- tbrca_data %>%
  filter(!is.na(concentration)) %>%
  mutate(pass_qc = ifelse(concentration >=dna_qc_threshold, "Yes", "No"))

samples_passing_qc <- nrow(tbrca_data_mod[tbrca_data_mod$pass_qc == "Yes",])

samples_failing_qc <- nrow(tbrca_data_mod[tbrca_data_mod$pass_qc == "No",])

fail_rate <- round((samples_failing_qc / (samples_passing_qc + samples_failing_qc)) * 100, 1)

ggplot(tbrca_data_mod, aes(x = disease, y = concentration)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  labs(y = "DNA concentration (ng/ul)",
       x = "Tumour BRCA referrals",
       title = paste0("DNA concentrations for ", nrow(tbrca_data_mod), " tumour BRCA referrals"),
       subtitle = paste0("Samples below ", dna_qc_threshold, " ng/ul: ", samples_failing_qc,
                         " (", fail_rate, "%)"),
       caption = paste0("Median DNA concentration: ", median(tbrca_data_mod$concentration), 
                        " ng/ul")) +
  ylim(0, 700) +
  geom_hline(yintercept = 3.3, linetype = "dashed")

##################################################
# Spread of Myriad GI Scores
##################################################

tbrca_data_collection <- read_excel(paste0(hrd_data_path, 
                                           "HRD TBRCA data collection Manchester_NEW_from Oct2022_2023.xlsx"),
                                    skip = 1) %>%
  janitor::clean_names()

tbrca_data_collection_clean <- tbrca_data_collection %>%
  filter(!gis_score_numerical_value_or_fail_not_tested %in% c("Fail", "Inconclusive",
                                                              "Not tested", NA)) %>%
  filter(t_brca_mutation_status != "Fail") %>%
  mutate(gi_score = as.numeric(gis_score_numerical_value_or_fail_not_tested),
         brca_mutation_clean = case_when(
           
           t_brca_mutation_status %in% c("Pathogenic BRCA1", "Pathogenic BRCA2") ~"BRCA positive",
           t_brca_mutation_status %in% c("No mutation detected", "no mutation detected") ~"BRCA negative"))

tbrca_data_collection_clean %>%
         group_by(gi_score, gis_pos_neg) %>%
         summarise(total = n()) %>%
  ggplot(aes(x = gi_score, y = total)) +
  #geom_col(aes(fill = gis_pos_neg), width = 1) +
  geom_smooth() +
  scale_fill_manual(values = c("#CCCCCC", "#FF6666")) +
  xlim(0, 100) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  #geom_vline(xintercept = 42, linetype = "dashed") +
  labs(title = paste0("Myriad GI score results for ", nrow(tbrca_data_collection_clean), 
                      " samples from the North West GLH"),
       y = "Number of samples",
       x = "Myriad GI score") 

compare_results %>%
  filter(!is.na(seqone_hrd_score)) %>%
  group_by(seqone_hrd_score) %>%
  summarise(total = n()) %>%
  ggplot(aes(x = seqone_hrd_score, y =total)) +
  geom_col(width = 0.01) +
  geom_smooth() +
  xlim(0, 1) +
  theme_bw() +
  theme(panel.grid = element_blank())

##################################################