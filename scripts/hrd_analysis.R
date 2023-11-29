# Homologous Recombination Deficiency Validation: Analysis
# joseph.shaw3@nhs.net; joseph.shaw2@mft.nhs.uk

# Setup -----------------------------------------------------------------------------

rm(list = ls())

## Packages -------------------------------------------------------------------------

library("ggpubr")
library("readxl")

## Functions ------------------------------------------------------------------------

source("functions/hrd_functions.R")

# Collate data ----------------------------------------------------------------------

## Collate Myriad data --------------------------------------------------------------

myriad_report_files <- list.files("data/myriad_reports/", full.names = TRUE)

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

seqone_report_files_v1_1 <- list.files("data/seqone_reports_v1_1/", 
                                  full.names = TRUE)

seqone_reports_v1_1 <- seqone_report_files_v1_1 |>
  map(\(seqone_report_files_v1_1) read_seqone_report(file = seqone_report_files_v1_1,
                                                version = "1.1")) |>
  list_rbind()

seqone_report_files_v1_2 <- list.files("data/seqone_reports_v1_2/", 
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
low_res_myriad_results <- read_excel(path =
  "data/myriad_reports_low_res/myriad_low_res_report_data.xlsx") |>
  mutate(nhs_number = as.numeric(gsub(pattern = "(\\D)", "", nhs_number)))

seraseq_gi_scores <- read_excel("data/seraseq_gi_scores.xlsx")

collated_myriad_info_mod <- rbind(
  collated_myriad_info, low_res_myriad_results,
  seraseq_gi_scores |>
    select(-dlms_dna_number)
)

## Samples DLMS information ---------------------------------------------------------

seqone_dlms_info <- read_csv(file = "data/seqone_dlms_info.csv")

## Pathology block ID check ---------------------------------------------------------

# Writing of pathology block IDs is inconsistent, so a manual check is required
# Seraseq controls classified as "pathology blocks match"

path_block_check <- read_csv(
  "data/pathology_block_id_check.csv",
  show_col_types = FALSE) 

## Initial DNA concentrations -------------------------------------------------------

# Exported from Sharepoint (HS2 Sample Prep 2023 - NEW.xlsx)
dna_concentrations <- read_excel(
  path = "data/dna_concentrations.xlsx",
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
# SSXT HS2 Library Prep 2023.xlsx

hs2_library_prep <- read_excel(
  path = "data/qpcr_qc.xlsx",
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

# Spreadsheet sent by Katie Sadler
# "HRD TBRCA data collection Manchester_NEW_from Oct2022_2023.xlsx"

tbrca_data_collection <- read_excel(
    "data/tbrca_data_collection.xlsx",
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

seqone_qc_data_v1.1 <- read_excel("data/seqone_qc_metrics_v1.1.xlsx") |>
  janitor::clean_names() |>
  dplyr::rename(
    shallow_sample_id = sample,
    read_length = read_len,
    insert_size = ins_size,
    million_reads = m_reads
  ) |> 
  mutate(version = "1.1")

seqone_qc_data_v1.2 <- read_excel("data/seqone_qc_metrics_v1.2.xlsx") |>
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
      select(shallow_sample_id, q_pcr_n_m, ts_ng_ul, total_yield, q_pcr_n_m,
             index, d1000_ts_size, ts_ng_ul, plate_position),
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

# Check all samples have a path block ID check and join information

samples_to_check <- join_tables |> 
  # Controls won't have all Myriad information - so remove before checking for NAs
  filter(!sample_type %in% c("Biobank control", "Seraseq control") &
           downsampled == "No") |> 
  # NCC column has NA values due 
  select(-ncc)

assert_that(anyNA.data.frame(samples_to_check) == FALSE)

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
                                                    "other")),
    
    robustness = as.numeric(robustness)
    
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

# Initial format - all DNA replicates included

v1.1_results <- compare_results |> 
  filter(path_block_manual_check == "pathology blocks match" & version == "1.1") |> 
  compare_tests()

v1.1_results_filtered <- compare_results |> 
  filter(path_block_manual_check == "pathology blocks match" & version == "1.1" &
           coverage >= coverage_threshold & input_ng >= input_threshold) |> 
  compare_tests()

v1.2_results <- compare_results |> 
  filter(path_block_manual_check == "pathology blocks match" & version == "1.2") |> 
  compare_tests()

v1.2_results_filtered <- compare_results |> 
  filter(path_block_manual_check == "pathology blocks match" & version == "1.2" &
           coverage >= coverage_threshold & input_ng >= input_threshold) |> 
  compare_tests()

metric_table_inputs <- rbind(add_version(v1.1_results[[2]], "SomaHRD v1.1"),
                      add_version(v1.1_results_filtered[[2]], "SomaHRD v1.1, thresholds applied"),
                      add_version(v1.2_results[[2]], "SomaHRD v1.2"),
                      add_version(v1.2_results_filtered[[2]], "SomaHRD v1.2, thresholds applied"))

export_timestamp(input = metric_table_inputs)

export_timestamp(input = v1.1_results[[1]])

export_timestamp(input = v1.2_results[[1]])

# New format - only 1 replicate per sample

v1.2_unique_sample_results <- compare_results |> 
  filter(path_block_manual_check == "pathology blocks match" & 
           version == "1.2") |> 
  filter(!duplicated(dlms_dna_number)) |> 
  compare_tests()

v1.2_unique_sample_results_filtered <- compare_results |> 
  filter(path_block_manual_check == "pathology blocks match" & 
           version == "1.2" &
           coverage >= coverage_threshold & 
           input_ng >= input_threshold) |> 
  filter(!duplicated(dlms_dna_number)) |> 
  compare_tests()

export_timestamp(input = v1.2_unique_sample_results[[1]])

metric_table_samples <- rbind(add_version(v1.2_unique_sample_results[[2]], 
                                          "SomaHRD v1.2"),
                        add_version(v1.2_unique_sample_results_filtered[[2]], 
                                    "SomaHRD v1.2, thresholds applied")) |> 
  select(-`DNA inputs`)

export_timestamp(input = metric_table_samples)

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
    y = "SeqOne HRD Probability (v1.2)") +
  ggpubr::stat_cor(method = "pearson", label.x = 60, label.y = 0.1) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 42, linetype = "dashed")

save_hrd_plot(input_height = 14, input_plot = hrd_score_plot)

## Robustness and low tumour content ------------------------------------------------

low_tumour_content_table <- compare_results |> 
  filter(version == "1.2" & low_tumour_fraction == "WARNING") |> 
  filter(!duplicated(dlms_dna_number)) |> 
  select(dlms_dna_number, low_tumour_fraction,
         seqone_hrd_status, myriad_hrd_status,
         path_block_manual_check, ncc, robustness) |> 
  arrange(dlms_dna_number) |> 
  rename("DNA Number" = dlms_dna_number,
         "Low tumour fraction status" = low_tumour_fraction, 
         "SeqOne HRD status (v1.2)" = seqone_hrd_status,
         "Myriad HRD status" = myriad_hrd_status,
         "Pathology block check" = path_block_manual_check,
         "NCC" = ncc,
         "Robustness" = robustness)
         
export_timestamp(input = low_tumour_content_table)

robustness_fail_samples <- compare_results |> 
  filter(version == "1.2" & robustness <= 0.85) 

robustness_warning_samples <- compare_results |> 
  filter(version == "1.2" & robustness < 0.93 &
           robustness > 0.85)

robustness_fail_table <- make_robustness_table(robustness_fail_samples)

robustness_warning_table <- make_robustness_table(robustness_warning_samples)

export_timestamp(input = robustness_fail_table)

export_timestamp(input = robustness_warning_table)

# Proportion of all samples with inconclusive results after input and coverage 
# filters applied
compare_results |> 
  filter(version == "1.2") |> 
  filter(!duplicated(dlms_dna_number)) |> 
  filter(coverage >= coverage_threshold & 
           input_ng >= input_threshold) |> 
  count(seqone_hrd_status) |> 
  mutate(percentage = round((n / sum(n))*100, 1))

## Intra-run variation --------------------------------------------------------------

intra_run_samples <- compare_results |> 
  filter(worksheet == "WS133557" & version == "1.2") |> 
  filter(base::duplicated(dlms_dna_number, fromLast = TRUE) |
           base::duplicated(dlms_dna_number, fromLast = FALSE)) 

intra_run_table <- intra_run_samples |> 
  select(sample_id, seqone_hrd_status, seqone_hrd_score, lga, lpc, coverage,
         robustness) |> 
  rename("Sample ID" = sample_id,
         "SeqOne HRD status (v1.2)" = seqone_hrd_status,
         "SeqOne HRD score (v1.2)"= seqone_hrd_score,
         "LGA" = lga,
         "LPC" = lpc,
         "Coverage (X)" = coverage,
         "Robustness" = robustness)

export_timestamp(input = intra_run_table)

intra_run_variation <- calculate_variation(intra_run_samples)

## Inter run variation --------------------------------------------------------------

inter_run_samples <- compare_results |>
  filter(version == "1.2" &
           # Remove intra-run replicates
           !sample_id %in% grep(x = compare_results$sample_id, pattern = "b|c",
                                value = TRUE)) |> 
  filter(base::duplicated(dlms_dna_number, fromLast = TRUE) |
    base::duplicated(dlms_dna_number, fromLast = FALSE))

lpc_lga_facet_plot <- plot_lpc_lga(inter_run_samples) +
  labs(x = "LGA", y = "LPC") +
  geom_point(data = inter_run_samples |>  
               filter(shallow_sample_id == "WS133557_21003549"),
             shape = 21, colour = safe_red, fill = NA, size = 4, stroke = 1) +
  facet_wrap(~dlms_dna_number)

save_hrd_plot(lpc_lga_facet_plot)

hrd_score_facet_plot <- ggplot(inter_run_samples, aes(x = worksheet, 
                                                       y = seqone_hrd_score)) +
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

hrd_score_reads_facet_plot <- ggplot(inter_run_samples, aes(x = million_reads,
                                                               y = seqone_hrd_score)) +
  geom_point(size = 2, alpha = 0.6,
             aes(colour = seqone_hrd_status)) +
  scale_colour_manual(name = "SeqOne HRD status",
                      values = c(safe_blue, safe_red, safe_grey)) +
  theme_bw() +
  xlim(0, 70) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  theme(legend.position = "bottom") +
  labs(x = "Millions of reads", y = "SeqOne HRD probability (v1.2)") +
  facet_wrap(~dlms_dna_number)

save_hrd_plot(hrd_score_reads_facet_plot)

## Intra-run variation boxplots -----------------------------------------------------

repeat_variation <- calculate_variation(inter_run_samples)
  
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
  ggplot(aes(x = myriad_gi_score, y = seqone_hrd_score)) +
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
         low_tumour_fraction, robustness) |> 
  rename("DNA number" = dlms_dna_number,
         "SeqOne HRD status (v1.2)" = seqone_hrd_status,
         "LGA" = lga,
         "LPC" = lpc,
         "Low tumour fraction status" = low_tumour_fraction,
         "Robustness" = robustness)

export_timestamp(input = biobank_control_results)

## QC metrics -----------------------------------------------------------------------

coverage_line <- geom_hline(yintercept = coverage_threshold, linetype = "dashed")

input_line <- geom_hline(yintercept = input_threshold, linetype = "dashed")

filtered_results <- compare_results |> 
  filter(path_block_manual_check == "pathology blocks match" & version == "1.2")

cov_v12 <- plot_qc(yvar = coverage, alpha_number = 1) +
  coverage_line +
  ylim(0, 2.5) +
  labs(x = "DNA input", y = "Coverage (X)")

input_v12 <- plot_qc(yvar = input_ng, alpha_number = 1) +
  ylim(0, 52) +
  input_line +
  labs(x = "DNA input", y = "DNA input (ng)")

coverage_input_plot <- ggarrange(cov_v12, input_v12,
          nrow = 2, ncol = 1,
          common.legend = TRUE,
          legend = "bottom")

save_hrd_plot(coverage_input_plot, input_width = 16)

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
         myriad_hrd_status, million_reads, coverage) |> 
  rename("DNA number" = dlms_dna_number,
         "Downsampled" = downsampled,
         "SeqOne HRD score (v1.2)" = seqone_hrd_score,
         "SeqOne HRD status (v1.2)" = seqone_hrd_status,
         "Myriad HRD status" = myriad_hrd_status,
         "Reads (millions)" = million_reads,
         "Coverage" = coverage)

export_timestamp(input = downsampled_table)

## Impact of SomaHRD version on HRD status ------------------------------------------

check_version_statuses <- compare_results |> 
  filter(downsampled == "No") |> 
  select(shallow_sample_id, seqone_hrd_status, version) |> 
  pivot_wider(id_cols = shallow_sample_id,
              names_from = version,
              values_from = seqone_hrd_status) |> 
  mutate(version_check = case_when(
    
    `1.1` == `1.2` ~"HRD statuses are the same",
    
    `1.1` != `1.2` ~"HRD statuses are NOT the same"
    
  ))

check_version_statuses |> 
  count(version_check)

# Manchester tBRCA service audit ----------------------------------------------------

## Manchester tBRCA DNA concentrations ----------------------------------------------

tbrca_data <- read.csv(file = "data/tbrca_dlms_info.csv")

dna_qc_threshold <- round(50 / 15, 0)

tbrca_data_mod <- tbrca_data |>
  filter(!is.na(concentration)) |>
  mutate(pass_qc = ifelse(concentration >= dna_qc_threshold, "Yes", "No"))

tbrca_data_mod |> 
  group_by(pass_qc) |> 
  summarise(total = n()) |> 
  mutate(proportion = round((total / sum(total))*100, 1))

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

tbrca_gi_scores |> 
  mutate(borderline = ifelse(gi_score <= 47 & gi_score >=37, "Yes", "No")) |> 
  count(borderline) |> 
  mutate(percentage = round((n / sum(n))*100, 1))

myriad_sd <- 2.2

borderline_myriad_samples <- tbrca_gi_scores |> 
  filter(gi_score <= 44 & gi_score >= 40)

## Manchester tBRCA inconclusive rate -----------------------------------------------

exclude <- c("Fail", "Not tested", "Awaiting result")

tbrca_inconclusive_data <- tbrca_data_collection_clean |> 
  filter(!gis_score_numerical_value_or_fail_not_tested %in% exclude &
           !gis_pos_neg %in% exclude &
           !overall_hrd_pos_neg %in% exclude) |> 
  mutate(result_type = case_when(
    gis_score_numerical_value_or_fail_not_tested == "Inconclusive" ~"Inconclusive",
    gi_score >=0 & gi_score <= 100 ~"Conclusive result"
  ))
  
tbrca_inconclusive_summary <- tbrca_inconclusive_data |> 
  count(result_type) |> 
  mutate(rate = round((n/sum(n))*100, 1))

# Export combined tables ------------------------------------------------------------

export_timestamp(input = compare_results)

# Inter-run variability investigation -----------------------------------------------

ggplot(inter_run_samples, aes(x = coverage,
                              y = lpc)) +
  geom_point(size = 2, alpha = 0.6,
             aes(colour = seqone_hrd_status)) +
  scale_colour_manual(name = "SeqOne HRD status",
                      values = c(safe_blue, safe_red, safe_grey)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  facet_wrap(~dlms_dna_number)

inter_run_mod <- inter_run_samples |> 
  mutate(total_bp = (million_reads*1000000) * read_length,
         coverage_calc = total_bp / (3.2 * 1000000000)) 

investigate_plot(lga, lpc)

investigate_plot(read_length, seqone_hrd_score)

investigate_plot(percent_aligned, percent_q30)

investigate_plot(index, million_reads)

investigate_plot(d1000_ts_size, lpc)

investigate_plot(percent_dups, million_reads)

investigate_plot(coverage_calc, coverage)

## Standard deviation ----------------------------------------------------------------

repeated_samples <- compare_results |>
  filter(version == "1.2") |> 
  filter(base::duplicated(dlms_dna_number, fromLast = TRUE) |
           base::duplicated(dlms_dna_number, fromLast = FALSE))

sd_score <- calculate_pooled_sd(repeated_samples, seqone_hrd_score)

sd_lga <- calculate_pooled_sd(repeated_samples, lga)

sd_lpc <- calculate_pooled_sd(repeated_samples, lpc)

sd_table <- tribble(
  ~Metric,     ~`Pooled standard deviation`, ~`Range of variation`,
  "HRD score", sd_score[[2]],                sd_score[[3]], 
  "LGA",       sd_lga[[2]],                  sd_lga[[3]],
  "LPC",       sd_lpc[[2]],                  sd_lpc[[3]]
)

export_timestamp(input = sd_table)

sd_results <- sd_score[[1]]

ggplot(sd_results, aes(x = reorder(dlms_dna_number, sd), y = sd)) +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "") +
  geom_hline(yintercept = sd_score[[2]], linetype = "dashed") +
  geom_hline(yintercept = 2* sd_score[[2]], linetype = "dashed")

# Remove outlier sample
calculate_pooled_sd(repeated_samples |>  filter(dlms_dna_number != "21003549"),
                    seqone_hrd_score,
                    round_places = 4)

calculate_pooled_sd(repeated_samples,
                    seqone_hrd_score,
                    round_places = 4)

hrd_score_repeat_plot <- repeated_samples |> 
  mutate(dlms_dna_number = factor(dlms_dna_number)) |> 
  mutate(dlms_dna_number = (fct_reorder(.f = dlms_dna_number, .x = seqone_hrd_score, 
                                        .fun = median))) |> 
  ggplot(aes(x = dlms_dna_number, y = seqone_hrd_score)) +
  geom_jitter(size = 4, alpha = 0.5, aes(colour = dlms_dna_number)) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "Sample", y = "SeqOne HRD Score")

repeated_samples |> 
  mutate(dlms_dna_number = factor(dlms_dna_number)) |> 
  mutate(dlms_dna_number = (fct_reorder(.f = dlms_dna_number, .x = seqone_hrd_score, 
                                        .fun = median))) |> 
  ggplot(aes(x = dlms_dna_number, y = seqone_hrd_score)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "Sample", y = "SeqOne HRD Score")

## Sample 21003549 ------------------------------------------------------------------

sample_21003549 <- read_excel("data/21003549_lga_lpc_counts.xlsx",
                              col_types = c("text", "numeric", "numeric", "numeric"))

ggplot(sample_21003549, aes(x = shallow_sample_id, y = lga, fill = shallow_sample_id)) +
  geom_col(position = "dodge") +
  theme_bw() +
  facet_wrap(~chromosome)

ggplot(sample_21003549, aes(x = shallow_sample_id, y = lpc, fill = shallow_sample_id)) +
  geom_col(position = "dodge") +
  theme_bw() +
  facet_wrap(~chromosome)

compare_results |> 
  filter(version == "1.2" & dlms_dna_number == 21003549) |> 
  select(shallow_sample_id, lga, lpc, seqone_hrd_score, robustness)

# Q: do samples with lower robustness tend to have lower LGA and LPC scores?

# A: no.

compare_results |> 
  filter(version == "1.2") |> 
  filter(duplicated(dlms_dna_number, fromLast = TRUE) |
           duplicated(dlms_dna_number, fromLast = FALSE)) |> 
  ggplot(aes(x = robustness, y = lga)) +
  geom_point() +
  facet_wrap(~dlms_dna_number)

compare_results |> 
  filter(version == "1.2") |> 
  filter(duplicated(dlms_dna_number, fromLast = TRUE) |
           duplicated(dlms_dna_number, fromLast = FALSE)) |> 
  ggplot(aes(x = robustness, y = lpc)) +
  geom_point() +
  facet_wrap(~dlms_dna_number)

## Genome profile check -------------------------------------------------------------

profile_check <- read_excel(path = 
                              "data/2023_11_20_08_22_11_genomic_profile_check.xlsx") |> 
  mutate(telomere_copy_profile = factor(telomere_copy_profile,
                                        levels = c("Normal", "Increased low",
                                                   "Increased high", "Decreased")))

telomere_colours = c("#CCCCCC", "#FF3366", "#CC0000", "#3300FF" )

telomere_profile_summary <- profile_check |> 
  count(telomere_copy_profile) 

export_timestamp(input = telomere_profile_summary)

results_and_profile <- compare_results |>
  filter(version == "1.2") |> 
  left_join(profile_check, by = "shallow_sample_id")

# Q: Are the strange telomere profiles seen on each worksheet?

# A: Decreased and increased high patterns are not seen on mid-output runs. But
# this may be due to smaller sample numbers.

ggplot(results_and_profile, aes(x = worksheet, 
                                y = )) +
  geom_bar(aes(fill = telomere_copy_profile), position = "dodge") +
  scale_fill_manual(values = telomere_colours)

# Q: are the strange telomere profiles caused by low coverage?

# A: there is no correlation with low coverage in repeated samples.

results_and_profile |> 
  filter(duplicated(dlms_dna_number, fromLast = TRUE) |
           duplicated(dlms_dna_number, fromLast = FALSE)) |> 
  ggplot(aes(x = worksheet, y = coverage)) +
  geom_jitter(aes(colour = telomere_copy_profile), size = 3, alpha = 0.6) +
  scale_colour_manual(values = telomere_colours) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "") +
  facet_wrap(~dlms_dna_number)

# A: there is no correlation with coverage in the full cohort.

ggplot(results_and_profile, aes(x = reorder(shallow_sample_id, 
                                            coverage), y = coverage)) +
  geom_jitter(aes(colour = telomere_copy_profile), size = 3, alpha = 0.6) +
  scale_colour_manual(values = telomere_colours) +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  labs(x = "")

# Q: are sample indexes involved?

# A: the same indexes have different telomere copy profiles, so it is unlikely
# to be indexes.

index_check <- results_and_profile |> 
  select(shallow_sample_id, index, telomere_copy_profile) |> 
  arrange(index)

# Q: are strange telomere copy profiles linked to:
#     - robustness
#     - read length
#     - input
#     - insert size
#     - percentage bases with quality over 30

# A: no, nothing obvious.

make_telomere_plot(robustness)

make_telomere_plot(read_length)

make_telomere_plot(input_ng)

make_telomere_plot(insert_size)

make_telomere_plot(percent_q30)

make_telomere_plot(percent_dups)

ggplot(results_and_profile, aes(x = percent_dups, y = percent_q30)) +
  geom_point(aes(colour = telomere_copy_profile)) +
  scale_colour_manual(values = telomere_colours) +
  theme_bw()

# Q: do we have a case where the same library had different telomere profiles?

# A: yes - sample 20112141

diff_profile <- results_and_profile |> 
  filter(dlms_dna_number == 20112141) |> 
  select(worksheet, sample_id, telomere_copy_profile,
         coverage, input_ng, lga, lpc)

# Q: if we ignore LPC and just look at LGA, is sample 21003549 still an outlier?

# A: yes.

results_and_profile |> 
  filter(duplicated(dlms_dna_number, fromLast = TRUE) |
  duplicated(dlms_dna_number, fromLast = FALSE)) |> 
  mutate(dlms_dna_number = as.character(dlms_dna_number)) |> 
  ggplot(aes(x = dlms_dna_number, y = lga)) +
  geom_point(size = 3, aes(colour = telomere_copy_profile), alpha = 0.5) +
  scale_colour_manual(values = telomere_colours) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "")

# Q: Were the samples that changed LGA and LPC results from v1.1 to v1.2 also the 
# ones with weird telomere results?

# A: No.

version_comparison <- compare_results |> 
  select(shallow_sample_id, dlms_dna_number, lga, lpc, version,
         seqone_hrd_score, seqone_hrd_status) |> 
  pivot_wider(id_cols = c(shallow_sample_id, dlms_dna_number),
              names_from = version,
              values_from = c(lga, lpc, seqone_hrd_score, seqone_hrd_status)) |> 
  mutate(lga_change = lga_1.2 - lga_1.1,
         lpc_change = lpc_1.2 - lpc_1.1,
         score_change = seqone_hrd_score_1.2 - seqone_hrd_score_1.1,
         status_change = ifelse(seqone_hrd_status_1.2 == seqone_hrd_status_1.1,
                                "Same", "Change")) |> 
  left_join(profile_check, by = "shallow_sample_id")
              
ggplot(version_comparison, aes(x = seqone_hrd_score_1.1, 
                    y = seqone_hrd_score_1.2)) +
  geom_jitter(aes(colour = telomere_copy_profile),
              size = 2, alpha = 0.6)

# Q: What were the inter-run results for v1.1 like?

# A: WS133557_21003549 was an outlier on v1.1 as well.

inter_run_1.1 <- compare_results |>
  filter(version == "1.1" &
           # Remove intra-run replicates
           !sample_id %in% grep(x = compare_results$sample_id, pattern = "b|c",
                                value = TRUE)) |> 
  filter(base::duplicated(dlms_dna_number, fromLast = TRUE) |
           base::duplicated(dlms_dna_number, fromLast = FALSE))

plot_lpc_lga(inter_run_1.1) +
  labs(x = "LGA", y = "LPC") +
  geom_point(data = inter_run_1.1 |>  
               filter(shallow_sample_id == "WS133557_21003549"),
             shape = 21, colour = safe_red, fill = NA, size = 4, stroke = 1) +
  facet_wrap(~dlms_dna_number)

# Q: what were the v1.1 results for 21003549?

# A: no change in LGA, LPC or score.

v1_1_check <- version_comparison |> 
  filter(dlms_dna_number == 21003549)

# Q: how many samples would adding a "telomere check" affect?

lga_frequency <- compare_results |> 
  filter(version == "1.2") |> 
  filter(!duplicated(dlms_dna_number)) |> 
  count(lga)

window <- lga_frequency |> 
  filter(lga >= 10 & lga <= 18)

lga_10_18_prop <- sum(window$n) / sum(lga_frequency$n)

# Proportion of cases with high or low telomere profiles
highlow_prop <- 10/98

# Predicted percentage of cases which would need repeating
(lga_10_18_prop * highlow_prop) * 100

# Q: Does robustness correlate with telomere profile for sample 21003549?

# A: no.

results_and_profile |> 
  filter(dlms_dna_number == 21003549) |> 
  ggplot(aes(x = robustness, y = telomere_copy_profile)) +
  geom_point(aes(colour = telomere_copy_profile),
              size = 2, alpha = 0.6)

# Q: was a sample with a high telomere copy profile ever run on a mid-output run?

# A: no.

results_and_profile |> 
  filter(telomere_copy_profile == "Increased high" & worksheet %in% c("WS135498",
                                                                      "WS134928"))

# Q: which samples had either increased or decreased telomere profiles?

# A: these ones

telomere_change <- results_and_profile |> 
  filter(telomere_copy_profile %in% c("Increased high", "Decreased")) |> 
  select(telomere_copy_profile,
         shallow_sample_id, lga, lpc, seqone_hrd_status, myriad_hrd_status,
         myriad_gi_score, 
         myriad_brca_status, myriad_patient_name,
         input_ng, coverage, path_block_manual_check) |>  
  arrange(telomere_copy_profile)

# Q: do we have any evidence of Myriad samples being repeat tested with different
# results?

# A: no. Only one sample with 2 Myriad results.

myriad_repeat <- tbrca_data_collection_clean |> 
  filter(duplicated(lab_number_id, fromLast = TRUE) |
           duplicated(lab_number_id, fromLast = FALSE)) |>  
  filter(!is.na(lab_number_id))

# Q: were the high telomere profile samples all in a similar plate location?

# A: no.

position_check <- results_and_profile |> 
  mutate(plate_y =sub(x = plate_position, pattern = "\\d", replacement = ""),
         plate_x = parse_number(plate_position)) |> 
  select(plate_position, telomere_copy_profile, shallow_sample_id) |> 
  arrange(plate_position)

# Are any of the biobank or seracare controls impacted?

control_check <- results_and_profile |> 
  filter(sample_type %in% c("Seraseq control", "Biobank control")) |> 
  select(shallow_sample_id, dlms_dna_number, 
         surname, telomere_copy_profile, seqone_hrd_status,
         lga, lpc) |> 
  arrange(dlms_dna_number)
