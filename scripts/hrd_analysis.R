# Homologous Recombination Deficiency Validation: Analysis
# joseph.shaw3@nhs.net

# Setup -----------------------------------------------------------------------------

rm(list = ls())

## Packages -------------------------------------------------------------------------

library("ggpubr")
library("readxl")

## Functions ------------------------------------------------------------------------

source("functions/hrd_functions.R")

# Collate data ----------------------------------------------------------------------

## Collate Myriad data --------------------------------------------------------------

myriad_report_files <- list.files(myriad_reports_location, full.names = TRUE)

collated_myriad_info <- myriad_report_files |>
  map(read_myriad_report) |>
  list_rbind()

stopifnot(setdiff(
  basename(myriad_report_files),
  collated_myriad_info$myriad_filename
) == 0)

# Change HRD status for BRCA positive sample with GIS of 5 (for SeqOne comparison later)

collated_myriad_info[collated_myriad_info$myriad_r_number == "R22-0LW4", 
                     "myriad_hrd_status"] <- "NEGATIVE"

## Collate SeqOne data --------------------------------------------------------------

seqone_report_files_v1_1 <- list.files(str_c(hrd_data_path, "seqone_reports_v1_1/"), 
                                  full.names = TRUE)

seqone_reports_v1_1 <- seqone_report_files_v1_1 |>
  map(\(seqone_report_files_v1_1) read_seqone_report(file = seqone_report_files_v1_1,
                                                version = "1.1")) |>
  list_rbind()

seqone_report_files_v1_2 <- list.files(str_c(hrd_data_path, "seqone_reports_v1_2/"), 
                                       full.names = TRUE)

seqone_reports_v1_2 <- seqone_report_files_v1_2 |>
  map(\(seqone_report_files_v1_2) read_seqone_report(file = seqone_report_files_v1_2,
                                                     version = "1.2")) |>
  list_rbind()

collated_seqone_info <- rbind(seqone_reports_v1_1, seqone_reports_v1_2)

# Check all files collated
stopifnot(setdiff(
  c(basename(seqone_report_files_v1_1), basename(seqone_report_files_v1_2)),
  collated_seqone_info$filename
) == 0)

# Check to make sure a report wasn't saved twice
stopifnot(nrow(collated_seqone_info |>
  # Combination of shallow sample ID and run date should be unique
  mutate(sample_date = str_c(shallow_sample_id, " ", date)) |>
  filter(duplicated(sample_date))) == 0)

## Edit Myriad information ----------------------------------------------------------

# pdf_text doesn't work on these reports as resolution is too low.
# I had to collate data manually in an Excel
low_res_myriad_results <- read_excel(path = paste0(
  hrd_data_path,
  "myriad_reports_low_res/myriad_low_res_report_data.xlsx"
)) |>
  mutate(nhs_number = as.numeric(gsub(pattern = "(\\D)", "", nhs_number)))

seraseq_gi_scores <- read_excel(paste0(hrd_data_path, "seraseq_gi_scores.xlsx"))

collated_myriad_info_mod <- rbind(
  collated_myriad_info, low_res_myriad_results,
  seraseq_gi_scores |>
    select(-dlms_dna_number)
)

## Samples DLMS information ---------------------------------------------------------

seqone_dlms_info <- read_csv(file = str_c(hrd_data_path, "seqone_dlms_info.csv"))

## Pathology block ID check ---------------------------------------------------------

# Writing of pathology block IDs is inconsistent, so a manual check is required
# Seraseq controls classified as "pathology blocks match"

path_block_check <- read_csv(
  paste0(
    hrd_data_path,
    "2023_10_24_14_10_20_export_for_check_edit.csv"
  ),
  show_col_types = FALSE) 

## Initial DNA concentrations -------------------------------------------------------

# Exported from Sharepoint
dna_concentrations <- read_excel(
  path = paste0(hrd_data_path, "HS2 Sample Prep 2023 - NEW.xlsx"),
  col_types = c(
    "date",    "numeric", "text",    "text",
    "date",    "numeric", "numeric", "numeric", 
    "numeric", "numeric", "date",    "numeric", 
    "numeric", "numeric", "numeric", "numeric",
    "date",    "text",    "numeric", "numeric", 
    "text",    "text",    "text",    "text"
  ),
  sheet = "HRD_SeqOne"
) |>
  janitor::clean_names() |> 
  dplyr::rename(
    dlms_dna_number = sample_id,
    forename = name,
    surname = x4,
    nanodrop_ng_ul = nanodrop,
    nd_260_280 = x7,
    nd_260_230 = x8,
    stock_volume_ul = qubit,
    qubit_dna_ul = x10,
    cleanup_date = clean_up,
    cleanup_nanodrop_ng_ul = x12,
    cleanup_260_280 = x13,
    cleanup_260_230 = x14,
    cleanup_volume = x15,
    cleanup_qubit_ng_ul = x16,
    dilution_date = dilution,
    dilution_worksheet = x18,
    dilution_dna_volume = x19,
    dilution_water_volume = x20
  )

dna_concentrations_mod <- dna_concentrations |>
  # Remove duplicated incomplete data
  slice(-(1:28)) |> 
  mutate(
    dilution_concentration = ifelse(

      comments %in% unique(grep(
        pattern = "neat",
        dna_concentrations$comments, value = TRUE
      )),
      qubit_dna_ul,
      round((dilution_dna_volume * qubit_dna_ul) /
        (dilution_water_volume + dilution_dna_volume), 2)
    ),

    # 15ul of diluted or neat (undiluted) DNA used in fragmentation reaction
    input_ng = dilution_concentration * 15,
    
    # Sequencing worksheet was not included on spreadsheet - add by date
    dilution_worksheet = case_when(
      
      as.character(date_submitted) %in% c("2023-07-21",
                                          "2023-07-26",
                                          "2023-08-03") ~"WS133557",
      
      as.character(date_submitted) == "2023-09-15" ~"WS134687",
      
      as.character(date_submitted) == "2023-09-27" ~"WS135001",
      
      as.character(date_submitted) == "2023-10-12" ~"WS135498")

  )

## qPCR library QC ------------------------------------------------------------------

# Copy of technologist spreadsheet originally saved at:
# S:/central shared/Genetics/Repository/Technical Teams/NGS/SureSelect XT HS2/

hs2_library_prep <- read_excel(
  path = str_c(hrd_data_path, "SSXT HS2 Library Prep 2023.xlsx"),
  sheet = "Sheet1",
  col_types = c(
    "text",    "guess",   "text",    "text",
    "text",    "numeric", "numeric", "numeric", 
    "numeric", "numeric", "numeric", "numeric", 
    "numeric", "numeric", "numeric", "numeric", 
    "numeric", "numeric", "text",    "guess", 
    "guess",   "guess",   "guess",   "guess", 
    "guess",   "guess",   "guess"
  )
) |>
  janitor::clean_names() |>
  filter(!plate_position %in% c("Plate Position", NA))

kapa_data_collated <- rbind(
  extract_kapa_data("133557", 20),
  extract_kapa_data("134687", 31),
  extract_kapa_data("134928", 7),
  extract_kapa_data("135001", 31),
  extract_kapa_data("135498", 7)
)

## tBRCA data -----------------------------------------------------------------------

tbrca_data_collection <- read_excel(
  paste0(
    hrd_data_path,
    "HRD TBRCA data collection Manchester_NEW_from Oct2022_2023.xlsx"
  ),
  skip = 1
) |>
  janitor::clean_names()

tbrca_data_collection_clean <- tbrca_data_collection |>
  mutate(
    gi_score = as.numeric(gis_score_numerical_value_or_fail_not_tested,
                          na = c(
                            "Fail", "Inconclusive",
                            "Not tested", NA
                          )),
    brca_mutation_clean = case_when(
      t_brca_mutation_status %in%
        c("Pathogenic BRCA1", "Pathogenic BRCA2") ~ "BRCA positive",
      t_brca_mutation_status %in%
        c("No mutation detected", "no mutation detected") ~ "BRCA negative"
    )
  )

## SeqOne QC data -------------------------------------------------------------------

seqone_qc_data_v1.1 <- read_excel(str_c(hrd_data_path, "seqone_qc_metrics_v1.1.xlsx")) |>
  janitor::clean_names() |>
  dplyr::rename(
    shallow_sample_id = sample,
    read_length = read_len,
    insert_size = ins_size,
    million_reads = m_reads
  ) |> 
  mutate(version = "1.1")

seqone_qc_data_v1.2 <- read_excel(str_c(hrd_data_path, "seqone_qc_metrics_v1.2.xlsx")) |>
  janitor::clean_names() |>
  dplyr::rename(
    shallow_sample_id = sample,
    read_length = read_len,
    insert_size = ins_size,
    million_reads = m_reads
  ) |> 
  mutate(version = "1.2")

seqone_qc_data <- rbind(seqone_qc_data_v1.1, seqone_qc_data_v1.2)

## Edit SeqOne information ----------------------------------------------------------

seqone_mod <- collated_seqone_info |>
  # Analysis on September 1st and before was previous pipeline version (1.0)
  filter(date > "2023-09-01") |>
  left_join(
    seqone_dlms_info |>
      select(
        dlms_dna_number, nhs_number, firstname, surname,
        i_gene_r_no, pathno, ncc
      ),
    by = "dlms_dna_number"
  ) |>
  mutate(
    downsampled = ifelse(grepl(pattern = c("downsampl|_0.5"), 
                         x = sample_id), "Yes", "No"),
    
    # WS134928 used the same libraries prepared on WS134687
    dilution_worksheet = ifelse(worksheet == "WS134928", "WS134687", worksheet)
    
  ) |>
  left_join(seqone_qc_data |> 
              select(-coverage), by = c("shallow_sample_id", "version"))

## Join data tables together --------------------------------------------------------

join_tables <- seqone_mod |>
  
  # Add Myriad results
  left_join(collated_myriad_info_mod, by = "nhs_number", keep = FALSE) |>
  
  # Add pathology block check
  left_join(path_block_check |> 
              select(dlms_dna_number, path_block_manual_check), 
            by = "dlms_dna_number",  keep = FALSE) |>
  
  # Add DNA concentrations
  left_join(dna_concentrations_mod |>
    select(dlms_dna_number, dilution_worksheet, qubit_dna_ul, input_ng), 
    by = c("dlms_dna_number", "dilution_worksheet"),
    keep = FALSE) |>
  
  # Add qPCR data
  left_join(
    kapa_data_collated |>
      select(shallow_sample_id, q_pcr_n_m, ts_ng_ul, total_yield, q_pcr_n_m),
    by = "shallow_sample_id",  keep = FALSE
  ) |> 
  
  # Add identities for validation table
  mutate(
    sample_type = case_when(
      
      grepl(pattern = "Biobank", x = surname, ignore.case = TRUE) ~ "Biobank control",
      
      grepl(pattern = "seraseq", x = surname, ignore.case = TRUE) ~ "Seraseq control",
      
      !grepl(pattern = "Biobank", x = surname, ignore.case = TRUE) &
        !grepl(pattern = "seraseq", x = surname, ignore.case = TRUE) &
        path_block_manual_check == path_block_match_text ~"Patient - same pathology block",
      
      !grepl(pattern = "Biobank", x = surname, ignore.case = TRUE) &
        !grepl(pattern = "seraseq", x = surname, ignore.case = TRUE) &
        path_block_manual_check == path_block_no_match_text ~"Patient - different pathology block",
      
      TRUE ~ "other"
      
    )
  )

# Analyse Data ----------------------------------------------------------------------

# Worksheets
worksheet_summary <- join_tables |> 
  filter(version == "1.2" & downsampled == "No") |> 
  group_by(worksheet) |> 
  count()

# DNA inputs
dna_input_summary <- join_tables |> 
  filter(version == "1.2" & downsampled == "No") |> 
  count(sample_type) |> 
  arrange(desc(n))

total_input_number <- sum(dna_input_summary$n)

Description <- c("Seraseq controls (Myriad GI scores available online)",
                       "MCRC biobank normal controls (no Myriad GI scores)",
                       "DNA samples from patients with DNA from a different pathology 
                       block tested by Myriad",
                       "DNA samples from patients with DNA from the same pathology 
                       block tested by Myriad")

# Samples
sample_summary <- join_tables |>
  filter(version == "1.2" & downsampled == "No") |>
  filter(!duplicated(dlms_dna_number)) |> 
  count(sample_type) |> 
  arrange(n) |> 
  cbind(Description) |> 
  select(Description, n) |> 
  rename("Number of samples" = n)

export_timestamp(input = sample_summary)

## Compare SeqOne and Myriad results ------------------------------------------------

compare_results <- join_tables |> 
  filter(downsampled == "No") |>
  mutate(
    
    outcome = case_when(
      
      seqone_hrd_status == pos_text & myriad_hrd_status == pos_text ~true_pos_text,
      
      seqone_hrd_status == neg_text & myriad_hrd_status == neg_text ~true_neg_text,
      
      seqone_hrd_status == pos_text & myriad_hrd_status == neg_text ~false_pos_text,
      
      seqone_hrd_status == neg_text & myriad_hrd_status == pos_text ~false_neg_text,
      
      seqone_hrd_status == incon_text & myriad_hrd_status == pos_text ~incon_pos_text,
      
      seqone_hrd_status == incon_text & myriad_hrd_status == neg_text ~incon_neg_text,
      
      TRUE ~ "other"),
    
    outcome_binary = case_when(
      
      outcome %in% c(true_pos_text, true_neg_text) ~consistent_text,
      
      outcome %in% c(false_pos_text, false_neg_text) ~inconsistent_text,
      
      outcome %in% c(incon_pos_text, incon_neg_text) ~inconclusive_text,
      
      TRUE ~ "other"),
    
    outcome_binary = fct(outcome_binary, levels = c(consistent_text,
                                                    inconsistent_text,
                                                    inconclusive_text,
                                                    "other"))
    
  )

# Every DNA input should have 2 rows: 1 for each pipeline version
assert_that(nrow(compare_results) == (total_input_number*2))

samples_without_myriad_results <- compare_results |> 
  filter(outcome == "other" | outcome_binary == "other")

# Only MCRC biobank controls should be without Myriad results
assert_that(unique(samples_without_myriad_results$sample_type) == "Biobank control")

## Sensitivity and specificity ------------------------------------------------------

coverage_threshold <- 0.5

input_threshold <- 47

v1.1_results <- compare_results |> 
  filter(path_block_manual_check == "pathology blocks match" & version == "1.1") |> 
  compare_tests(outcome)

v1.1_results_filtered <- compare_results |> 
  filter(path_block_manual_check == "pathology blocks match" & version == "1.1" &
           coverage >= coverage_threshold & input_ng >= input_threshold) |> 
  compare_tests(outcome)

v1.2_results <- compare_results |> 
  filter(path_block_manual_check == "pathology blocks match" & version == "1.2") |> 
  compare_tests(outcome)

v1.2_results_filtered <- compare_results |> 
  filter(path_block_manual_check == "pathology blocks match" & version == "1.2" &
           coverage >= coverage_threshold & input_ng >= input_threshold) |> 
  compare_tests(outcome)

metric_table <- rbind(add_version(v1.1_results[[2]], "SomaHRD v1.1"),
                      add_version(v1.1_results_filtered[[2]], "SomaHRD v1.1, thresholds applied"),
                      add_version(v1.2_results[[2]], "SomaHRD v1.2"),
                      add_version(v1.2_results_filtered[[2]], "SomaHRD v1.2, thresholds applied"))

export_timestamp(hrd_output_path, metric_table)

export_timestamp(hrd_output_path, v1.1_results[[1]])

export_timestamp(hrd_output_path, v1.2_results[[1]])

## Myriad and SeqOne score correlation ----------------------------------------------

hrd_score_plot <- compare_results |> 
  filter(path_block_manual_check == "pathology blocks match" & 
           version == "1.2" &
           coverage >= coverage_threshold & 
           input_ng >= input_threshold) |> 
  ggplot(aes(myriad_gi_score, seqone_hrd_score)) +
  geom_point(size = 3, alpha = 0.7,
             aes(colour = seqone_hrd_status)) +
  scale_colour_manual(name = "SeqOne HRD Status",
                      values = c(safe_blue, safe_red, safe_grey)) +
  theme_bw() +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = c(0, 25, 42, 50, 75, 100)
  ) +
  ylim(0, 1) +
  theme(panel.grid = element_blank(),
        legend.position = "bottom") +
  labs(
    x = "Myriad Genome Instability Score",
    y = "SeqOne HRD Score") +
  ggpubr::stat_cor(method = "pearson", label.x = 60, label.y = 0.1) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 42, linetype = "dashed")

save_hrd_plot(input_height = 14, input_plot = hrd_score_plot)

## Robustness and low tumour content ------------------------------------------------

ltc_table <- compare_results |> 
  filter(version == "1.2" & low_tumour_fraction == "WARNING") |> 
  filter(!duplicated(dlms_dna_number)) |> 
  select(dlms_dna_number, low_tumour_fraction, seqone_hrd_status,
         myriad_hrd_status,
         path_block_manual_check, ncc, robustness) |> 
  arrange(dlms_dna_number)

export_timestamp(input = ltc_table)

robustness_fail_samples <- compare_results |> 
  filter(version == "1.2" & robustness <= 0.85) 

robustness_warning_samples <- compare_results |> 
  filter(version == "1.2" & robustness < 0.93 &
           robustness > 0.85)

robustness_fail_table <- make_robustness_table(robustness_fail_samples)

robustness_warning_table <- make_robustness_table(robustness_warning_samples)

export_timestamp(input = robustness_fail_table)

export_timestamp(input = robustness_warning_table)

## Intra-run variation --------------------------------------------------------------

intra_run_table <- compare_results |> 
  filter(worksheet == "WS133557" & version == "1.2") |> 
  filter(base::duplicated(dlms_dna_number, fromLast = TRUE) |
           base::duplicated(dlms_dna_number, fromLast = FALSE)) |> 
  select(sample_id, seqone_hrd_status, seqone_hrd_score, lga, lpc, coverage,
         robustness)

export_timestamp(input = intra_run_table)

## Inter run variation --------------------------------------------------------------

repeat_results_1_2 <- compare_results |>
  filter(version == "1.2" &
           # Remove intra-run replicates
           !sample_id %in% grep(x = compare_results$sample_id, pattern = "b|c",
                                value = TRUE)) |> 
  filter(base::duplicated(dlms_dna_number, fromLast = TRUE) |
    base::duplicated(dlms_dna_number, fromLast = FALSE))

lpc_lga_facet_plot <- plot_lpc_lga(repeat_results_1_2) +
  labs(x = "LGA", y = "LPC") +
  facet_wrap(~dlms_dna_number)

save_hrd_plot(lpc_lga_facet_plot)

hrd_score_facet_plot <- ggplot(repeat_results_1_2, aes(x = worksheet, y = seqone_hrd_score)) +
  geom_point(size = 2, 
             aes(colour = seqone_hrd_status)) +
  scale_colour_manual(name = "SeqOne HRD status",
                      values = c(safe_blue, safe_red, safe_grey)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "bottom") +
  labs(x = "", y = "SeqOne HRD probability") +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  facet_wrap(~dlms_dna_number) 

save_hrd_plot(hrd_score_facet_plot)

repeat_variation <- repeat_results_1_2 |> 
  group_by(dlms_dna_number) |> 
  summarise(max_lga = max(lga),
            min_lga = min(lga),
            range_lga = max_lga-min_lga,
            max_lpc = max(lpc),
            min_lpc = min(lpc),
            range_lpc = max_lpc - min_lpc,
            max_cov = max(coverage),
            min_cov = min(coverage),
            range_cov = max_cov - min_cov,
            max_score = max(seqone_hrd_score),
            min_score = min(seqone_hrd_score),
            range_score = max_score - min_score)

lga_var <- plot_variation(yvar = range_lga) +
  labs(y = "Difference in LGA") +
  ylim(0,30)

lpc_var <- plot_variation(yvar = range_lpc) +
  labs(y = "Difference in LPC") +
  ylim(0,30)

score_var <- plot_variation(yvar = range_score) +
  labs(y = "Difference in HRD score") +
  ylim(0, 1)

lga_lpc_variation <- ggarrange(lga_var, lpc_var, score_var, 
                               ncol = 3, nrow = 1)

save_hrd_plot(lga_lpc_variation, input_height = 6)

## Seraseq controls -----------------------------------------------------------------

seraseq_control_data <- compare_results |>
  filter(sample_type == "Seraseq control") |>
  mutate(firstname_factor = factor(firstname, levels = c(
    "FFPE HRD Negative",
    "Low-Positive FFPE HRD",
    "High-Positive FFPE HRD"
  )))

seraseq_plot <- seraseq_control_data |> 
  filter(version == "1.2") |> 
  ggplot(aes(x = myriad_gi_score, 
                                 y = seqone_hrd_score)) +
  geom_point(aes(colour = seqone_hrd_status),
             size = 2, alpha = 0.6) +
  scale_colour_manual(name = "SeqOne HRD Status",
                      values = c(safe_blue, safe_red, safe_grey)) +
  
  theme_bw() +
  theme(legend.position = "bottom") +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 42, linetype = "dashed") +
  ylim(0, 1) +
  xlim(0, 80) +
  labs(x = "SeqOne HRD score", y = "Myriad GI score") +
  facet_wrap(~firstname_factor)

save_hrd_plot(seraseq_plot, input_height = 10)

## Biobank controls -----------------------------------------------------------------

biobank_control_results <- compare_results |> 
  filter(sample_type == "Biobank control" & version == "1.2") |> 
  select(dlms_dna_number, seqone_hrd_status, lga, lpc,
         low_tumor_fraction, robustness)

export_timestamp(input = biobank_control_results)

## QC metrics -----------------------------------------------------------------------

coverage_line <- geom_hline(yintercept = coverage_threshold, linetype = "dashed")

input_line <- geom_hline(yintercept = input_threshold, linetype = "dashed")

filtered_results <- compare_results |> 
  filter(path_block_manual_check == "pathology blocks match" & version == "1.2") |> 
  mutate(robustness = as.numeric(robustness))

cov_v12 <- plot_qc(df = filtered_results, yvar = coverage, outcome = outcome_binary) +
  coverage_line +
  ylim(0, 2.5) +
  labs(x = "DNA input", y = "Coverage (X)")

input_v12 <- plot_qc(df = filtered_results, yvar = input_ng, outcome = outcome_binary) +
  ylim(0, 52) +
  input_line +
  labs(x = "DNA input", y = "DNA input (ng)")

coverage_input_plot <- ggarrange(cov_v12, input_v12,
          nrow = 2, ncol = 1,
          common.legend = TRUE,
          legend = "bottom")

save_hrd_plot(coverage_input_plot)

readlength_plot_v2 <- plot_qc(yvar = read_length) +
  ylim(50, 150) +
  labs(y = "Read length (bp)")

millreads_plot_v2 <- plot_qc(yvar = million_reads) +
  ylim(0, 80) +
  labs(y = "Millions of reads")

insertsize_plot_v2 <- plot_qc(yvar = insert_size) +
  ylim(0, 200) +
  labs(y = "Insert size (bp)")

percentq30_plot_v2 <- plot_qc(yvar = percent_q30) +
  ylim(80, 100) +
  labs(y = "Percentage reads Q30 (%)")

percentaligned_plot_v2 <- plot_qc(yvar = percent_aligned) +
  ylim(99, 101) +
  labs(y = "Percentage reads aligned (%)")

percentdups_plot_v2 <- plot_qc(yvar = percent_dups) +
  ylim(0, 8) +
  labs(y = "Percentage duplicate reads (%)")

total_yield_plot_v2 <- plot_qc(yvar = total_yield) +
  ylim(0, 2500) +
  labs(y = "qPCR total yield (ng)")

robustness_plot_v2 <- plot_qc(yvar = robustness) +
  ylim(0.5, 1) +
  labs(y = "Robustness")

qc_metrics_v2 <- ggarrange(readlength_plot_v2, millreads_plot_v2, 
                           insertsize_plot_v2, percentq30_plot_v2, 
                           percentaligned_plot_v2, percentdups_plot_v2,
                           total_yield_plot_v2, robustness_plot_v2,
                           ncol = 2, nrow = 4,
                           common.legend = TRUE,
                           legend = "bottom")

save_hrd_plot(qc_metrics_v2, input_height = 24, input_width = 16)

## LGA and LPC ----------------------------------------------------------------------

lga_vs_lpc <- compare_results |> 
  filter(version == "1.2") |> 
  plot_lpc_lga() +
  labs(
    x = "Large Genomic Alterations",
    y = "Loss of Parental Copy")

save_hrd_plot(lga_vs_lpc)


## Downsampling ---------------------------------------------------------------------

downsampled_samples <- grep(pattern = "20103853b|20112141b",
                      x = join_tables$sample_id,
                      value = TRUE)


downsampled_table <- join_tables |> 
  filter(sample_id %in% downsampled_samples & version == "1.2") |> 
  select(dlms_dna_number, downsampled, seqone_hrd_score, seqone_hrd_status, 
         myriad_hrd_status, million_reads, coverage)

export_timestamp(input = downsampled_table)

# Manchester tBRCA service audit ----------------------------------------------------

## Manchester tBRCA DNA concentrations ----------------------------------------------

tbrca_data <- read.csv(file = str_c(hrd_data_path, "tBRCA_dlms_info.csv"))

dna_qc_threshold <- round(50 / 15, 0)

tbrca_data_mod <- tbrca_data |>
  filter(!is.na(concentration)) |>
  mutate(pass_qc = ifelse(concentration >= dna_qc_threshold, "Yes", "No"))

tbrca_data_mod |> 
  group_by(pass_qc) |> 
  summarise(total = n(),
    proportion = total / sum(total))

samples_passing_qc <- nrow(tbrca_data_mod[tbrca_data_mod$pass_qc == "Yes", ])

samples_failing_qc <- nrow(tbrca_data_mod[tbrca_data_mod$pass_qc == "No", ])

fail_rate <- round((samples_failing_qc /
  (samples_passing_qc + samples_failing_qc)) * 100, 1)

## Manchester tBRCA GI scores -------------------------------------------------------

tbrca_gi_scores <- tbrca_data_collection_clean |>
  filter(!is.na(gi_score))

myriad_gi_profile_plot <- tbrca_gi_scores |> 
  ggplot(aes(x = gi_score, y = )) +
  geom_histogram(binwidth = 1, aes(fill = gis_pos_neg)) +
  scale_fill_manual(values = c(safe_blue, safe_red)) +
  theme_bw() +
  scale_x_continuous(
    breaks = c(0, 25, 42, 50, 75, 100)
  ) +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom"
  ) +
  labs(
    y = "Number of samples",
    x = "Myriad GI score",
    fill = "Myriad GI Status")

save_hrd_plot(myriad_gi_profile_plot)

## Manchester tBRCA inconclusive rate -----------------------------------------------

exclude <- c("Fail", "Not tested", "Awaiting result")

tbrca_inconclusive_data <- tbrca_data_collection_clean |> 
  filter(!gis_score_numerical_value_or_fail_not_tested %in% exclude &
           !gis_pos_neg %in% exclude &
           !overall_hrd_pos_neg %in% exclude) |> 
  mutate(result_type = case_when(
    gis_score_numerical_value_or_fail_not_tested == "Inconclusive" ~"Inconclusive",
    gi_score >=0 & gi_score <= 100 ~"Conclusive result"
  )) |> 
  count(result_type)

# Export combined tables ------------------------------------------------------------

export_timestamp(input = compare_results)
