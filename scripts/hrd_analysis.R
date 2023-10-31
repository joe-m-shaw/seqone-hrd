# Homologous Recombination Deficiency Validation: Analysis
# joseph.shaw3@nhs.net

# Setup -----------------------------------------------------------------------------

rm(list = ls())

## Packages and filepaths -----------------------------------------------------------

library("ggpubr")
library("readxl")
library("epiR")

## DLMS connection and functions ----------------------------------------------------

source("scripts/dlms_connection.R")
source("functions/hrd_functions.R")

# Collate data ----------------------------------------------------------------------

## Source scripts and collate data --------------------------------------------------

myriad_report_files <- list.files(myriad_reports_location, full.names = TRUE)

collated_myriad_info <- myriad_report_files |> 
  map(read_myriad_report) |> 
  list_rbind()

# Check all files included
stopifnot(setdiff(basename(myriad_report_files), collated_myriad_info$myriad_filename) == 0)

collated_seqone_info <- collate_seqone_reports()

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

## Pathology block ID check ---------------------------------------------------------

# Extract data from DLMS via ODBC
seqone_dlms_info <- get_sample_data(collated_seqone_info$dlms_dna_number) |>
  rename(
    dlms_dna_number = labno,
    nhs_number = nhsno
  ) |>
  mutate(
    ncc_single = sub(
      x = comments,
      pattern = ".+(\\W{1}\\d{2}%).+",
      replacement = "\\1"
    ),
    ncc_range = sub(
      x = comments,
      pattern = ".+(\\d{2}\\W{1}\\d{2}%).+",
      replacement = "\\1"
    ),
    ncc = ifelse(str_length(ncc_single) == 4 & ncc_range > 6,
      ncc_single, ifelse(str_length(ncc_range) == 6,
        ncc_range, NA
      )
    )
  )

stopifnot(setdiff(
  seqone_dlms_info$dlms_dna_number,
  collated_seqone_info$dlms_dna_number) == 0)

# Enter a fake NHS number for the Seraseq and Biobank controls, to allow
# Myriad scores to be joined later, and to keep Biobank controls within the
# dataset

seqone_dlms_info[seqone_dlms_info$dlms_dna_number == 23032086, "nhs_number"] <- 1
seqone_dlms_info[seqone_dlms_info$dlms_dna_number == 23032088, "nhs_number"] <- 2
seqone_dlms_info[seqone_dlms_info$dlms_dna_number == 23031639, "nhs_number"] <- 3
seqone_dlms_info[seqone_dlms_info$dlms_dna_number == 23033285, "nhs_number"] <- 4
seqone_dlms_info[seqone_dlms_info$dlms_dna_number == 23033279, "nhs_number"] <- 5
seqone_dlms_info[seqone_dlms_info$dlms_dna_number == 23033288, "nhs_number"] <- 6

export_for_check <- seqone_dlms_info |>
  filter(dlms_dna_number %in% collated_seqone_info$dlms_dna_number) |>
  left_join(collated_myriad_info_mod |>
    select(
      nhs_number, myriad_pathology_block_pg1,
      myriad_pathology_block_pg2
    ), by = "nhs_number") |>
  select(
    dlms_dna_number, firstname, surname, nhs_number,
    pathno, myriad_pathology_block_pg1, myriad_pathology_block_pg2
  ) |>
  arrange(nhs_number)

export_timestamp(hrd_data_path, export_for_check)

# Writing of pathology block IDs is inconsistent, so a manual check is required
# Seraseq controls classified as "pathology blocks match"

path_block_check <- read_csv(
  paste0(
    hrd_data_path,
    "2023_10_24_14_10_20_export_for_check_edit.csv"
  ),
  show_col_types = FALSE
) |>
  select(dlms_dna_number, path_block_manual_check)

## Initial DNA concentrations -------------------------------------------------------

dna_concentrations <- read_excel(
  path = paste0(hrd_data_path, "HS2_sample_prep_export.xlsx"),
  col_types = c(
    "date", "numeric", "text", "text",
    "date", "numeric", "numeric",
    "numeric", "numeric", "numeric",
    "date", "text", "numeric",
    "numeric", "text", "text"
  )
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
    dilution_date = dilution,
    dilution_worksheet = x12,
    dilution_dna_volume = x13,
    dilution_water_volume = x14
  ) |>
  filter(!is.na(dlms_dna_number))

dna_concentrations_mod <- dna_concentrations |>
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

    # 1 haploid genome = 3.3 picograms
    input_genomes = (input_ng * 1000) / 3.3
  ) |>
  # Some samples used on previous runs. Filtering by comments isolates
  # the dilutions for the HRD runs
  filter(!is.na(comments)) |>
  # Checked: repeat rows of the same DNA number have the same DNA concentration
  # and dilution volumes
  filter(!base::duplicated(dlms_dna_number))

## qPCR library QC ------------------------------------------------------------------

tech_team <- "S:/central shared/Genetics/Repository/Technical Teams/"

hs2_library_prep <- read_excel(
  path = paste0(tech_team, "NGS/SureSelect XT HS2/SSXT HS2 Library Prep 2023.xlsx"),
  sheet = "Sheet1",
  col_types = c(
    "text", "guess", "text", "text",
    "text", "numeric", "numeric",
    "numeric", "numeric", "numeric",
    "numeric", "numeric", "numeric",
    "numeric", "numeric", "numeric",
    "numeric", "numeric", "text",
    "guess", "guess", "guess", "guess",
    "guess", "guess", "guess", "guess"
  )
) |>
  janitor::clean_names() |>
  filter(!plate_position %in% c("Plate Position", NA))

kapa_data_collated <- rbind(
  extract_kapa_data("133557", 20),
  extract_kapa_data("134687", 31),
  extract_kapa_data("134928", 7),
  extract_kapa_data("135001", 31)
)

## SeqOne QC data -------------------------------------------------------------------

seqone_qc_data <- read_excel(paste0(
  hrd_data_path,
  "seqone_qc_metrics/seqone_qc_metrics_2023_10_03.xlsx"
)) |>
  janitor::clean_names() |>
  dplyr::rename(
    shallow_sample_id = sample,
    read_length = read_len,
    insert_size = ins_size,
    million_reads = m_reads
  )

## Amended SeqOne scores ------------------------------------------------------------

seqone_amended <- read_excel(
  path =
    paste0(
      hrd_data_path,
      "SeqOne_Manchester Amended Results 241023.xlsx"
    )
) |>
  janitor::clean_names() |>
  rename(
    shallow_sample_id = sample,
    seqone_hrd_status_amended = amended_seq_one_hrd_status,
    lga_amended = lga,
    lpc_amended = lpc
  )


## Edit SeqOne information ----------------------------------------------------------

seqone_mod <- collated_seqone_info |>
  # Analysis on September 1st and before was previous pipeline version

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
    downsampled = ifelse(grepl(pattern = "downsampl", x = sample_id), "Yes", "No"),
    identity = ifelse(grepl(
      pattern = "Biobank|seraseq", x = surname,
      ignore.case = TRUE
    ),
    "Control", "Patient"
    ),
    control_type = case_when(
      grepl(pattern = "Biobank", x = surname, ignore.case = TRUE) ~ "Biobank control",
      grepl(pattern = "seraseq", x = surname, ignore.case = TRUE) ~ "Seraseq control",
      TRUE ~ "patient"
    )
  ) |>
  left_join(seqone_qc_data, by = "shallow_sample_id")

## Join data tables together --------------------------------------------------------

consistent_text <- "Seqone HRD status consistent with Myriad"

inconsistent_text <- "Seqone HRD status NOT consistent with Myriad"

inconclusive_text <- "SeqOne HRD status inconclusive"

compare_results <- seqone_mod |>
  left_join(collated_myriad_info_mod, by = "nhs_number") |>
  left_join(path_block_check, by = "dlms_dna_number") |>
  mutate(hrd_status_check = case_when(
    myriad_hrd_status == seqone_hrd_status ~ consistent_text,
    myriad_hrd_status != seqone_hrd_status ~ inconsistent_text,
    TRUE ~ "other"
  )) |>
  filter(downsampled == "No") |>
  left_join(dna_concentrations_mod |>
    select(
      dlms_dna_number, qubit_dna_ul, input_ng,
      input_genomes
    ), by = "dlms_dna_number") |>
  left_join(
    kapa_data_collated |>
      select(shallow_sample_id, q_pcr_n_m, ts_ng_ul, total_yield, q_pcr_n_m),
    by = "shallow_sample_id"
  ) |>
  left_join(
    seqone_amended |>
      select(
        shallow_sample_id, seqone_hrd_status_amended,
        lga_amended, lpc_amended, low_tumor_fraction_returns_a_warning_only
      ),
    by = "shallow_sample_id"
  ) |>
  mutate(
    change_in_status = case_when(
      seqone_hrd_status == seqone_hrd_status_amended ~ "No change",
      seqone_hrd_status != seqone_hrd_status_amended ~ "Change"
    ),
    hrd_status_check_amended = case_when(
      myriad_hrd_status == seqone_hrd_status_amended &
        seqone_hrd_status_amended != "NON-CONCLUSIVE" ~ consistent_text,
      myriad_hrd_status != seqone_hrd_status_amended &
        seqone_hrd_status_amended != "NON-CONCLUSIVE" ~ inconsistent_text,
      myriad_hrd_status != seqone_hrd_status_amended &
        seqone_hrd_status_amended == "NON-CONCLUSIVE" ~ inconclusive_text,
      TRUE ~ "other"
    )
  )

# Analyse Data ----------------------------------------------------------------------

## Inter run variation --------------------------------------------------------------

repeat_results <- compare_results |>
  filter(base::duplicated(dlms_dna_number, fromLast = TRUE) |
    base::duplicated(dlms_dna_number, fromLast = FALSE)) |>
  filter(downsampled == "No") |>
  mutate(input_category = case_when(
    input_ng >= 48 ~ "50ng input",
    input_ng < 48 ~ "lower than 50ng input"
  ))

repeat_facet_plot <- ggplot(repeat_results, aes(
  x = worksheet,
  y = seqone_hrd_score
)) +
  geom_point(
    alpha = 0.5,
    aes(
      shape = seqone_hrd_status,
      size = coverage.x
    )
  ) +
  facet_wrap(~dlms_dna_number) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90),
    legend.position = "bottom"
  ) +
  labs(
    title = "SeqOne results for repeated samples",
    x = "",
    y = "SeqOne HRD score",
    caption = "Plot name: repeat_facet_plot"
  ) +
  geom_hline(yintercept = 0.50, linetype = "dashed")

save_hrd_plot(repeat_facet_plot, input_width = 15, input_height = 15)

sample_20127786_plot <- make_individual_plot(20127786)

sample_21011999_plot <- make_individual_plot(21011999)

sample_23032088_plot <- make_individual_plot(23032088)

sample_21013520_plot <- make_individual_plot(21013520)

## Intra-run variation --------------------------------------------------------------

intra_run_plot <- compare_results |>
  filter(worksheet == "WS133557") |>
  filter(base::duplicated(dlms_dna_number, fromLast = TRUE) |
    base::duplicated(dlms_dna_number, fromLast = FALSE)) |>
  ggplot(aes(
    x = sample_id,
    y = seqone_hrd_score
  )) +
  geom_point(size = 3) +
  labs(
    x = "", y = "SeqOne HRD Score",
    title = "Repeated Sample Results",
    subtitle = "Worksheet WS133557"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))

save_hrd_plot(intra_run_plot, input_width = 10, input_height = 10)

## Repeat testing summaries ---------------------------------------------------------

summary_21013520 <- get_sample_summary_info(21013520)

export_timestamp(hrd_output_path, summary_21013520)

summary_23032088 <- get_sample_summary_info(23032088)

export_timestamp(hrd_output_path, summary_23032088)

summary_20127786 <- get_sample_summary_info(20127786)

export_timestamp(hrd_output_path, summary_20127786)

summary_21011999 <- get_sample_summary_info(21011999)

export_timestamp(hrd_output_path, summary_21011999)

## Compare results of mid output run - WS134928 -------------------------------------

mid_output_run <- seqone_mod |>
  filter(worksheet == "WS134928")

WS134928_repeated_samples <- mid_output_run$dlms_dna_number

WS134928_WS134687_data <- compare_results |>
  filter(dlms_dna_number %in% WS134928_repeated_samples &
    worksheet %in% c("WS134687", "WS134928"))

WS134928_WS134687_plot <- ggplot(
  WS134928_WS134687_data,
  aes(x = coverage.x, y = seqone_hrd_score)
) +
  geom_point(size = 3, aes(
    shape = seqone_hrd_status,
    colour = worksheet
  )) +
  facet_wrap(~dlms_dna_number) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(
    x = "Coverage", y = "SeqOne HRD Score",
    title = "Sample libraries repeated on mid-output run (WS134928)",
    caption = "Plot name: WS134928_WS134687_plot"
  ) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  xlim(0.8, 1.8)

ggplot(
  WS134928_WS134687_data,
  aes(x = lga, y = coverage.x)
) +
  geom_point(size = 3, alpha = 0.5, aes(
    shape = seqone_hrd_status,
    colour = worksheet
  )) +
  facet_wrap(~dlms_dna_number) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(
    title = "Sample libraries repeated on mid-output run (WS134928)",
    caption = "Plot name: WS134928_WS134687_plot"
  )

## HRD status comparison with Myriad ------------------------------------------------

comparison_summary <- compare_results |>
  filter(!is.na(myriad_gi_score) & myriad_gi_score != "NA") |>
  filter(path_block_manual_check != "NA") |>
  group_by(hrd_status_check, path_block_manual_check) |>
  summarise(total = n())

myriad_comparison_plot <- ggplot(comparison_summary, aes(
  x = hrd_status_check,
  y = total
)) +
  geom_col() +
  geom_text(aes(label = total), vjust = -0.5) +
  theme_bw() +
  labs(
    y = "Total DNA Inputs", x = "",
    title = "Consistency of Myriad and Seqone Testing",
    subtitle = paste0("Data for ", sum(comparison_summary$total), " DNA inputs"),
    caption = "Plot name: myriad_comparison_plot"
  ) +
  facet_wrap(~path_block_manual_check) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Impact of pathology blocks -------------------------------------------------------

results_for_path_block_plot <- compare_results |>
  filter(path_block_manual_check != "NA" & !is.na(myriad_gi_score))

path_block_plot <- ggplot(
  results_for_path_block_plot,
  aes(x = myriad_gi_score, y = seqone_hrd_score)
) +
  geom_point(size = 3, alpha = 0.6, aes(shape = hrd_status_check)) +
  theme_bw() +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = c(0, 25, 42, 50, 75, 100)
  ) +
  ylim(0, 1) +
  theme(panel.grid = element_blank()) +
  labs(
    x = "Myriad Genome Instability Score",
    y = "SeqOne HRD Score",
    title = "Comparison of Myriad vs SeqOne HRD Testing",
    subtitle = paste0("Data for ", nrow(results_for_path_block_plot), " DNA inputs")
  ) +
  ggpubr::stat_cor(method = "pearson", label.x = 50, label.y = 0.25) +
  # facet_wrap(~path_block_manual_check) +
  theme_bw() +
  theme(
    panel.grid = element_blank(), legend.position = "bottom",
    legend.title = element_blank()
  ) +
  guides(shape = guide_legend(ncol = 1)) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 42, linetype = "dashed")

high_quality_results <- compare_results |>
  filter(path_block_manual_check != "NA" & !is.na(myriad_gi_score)) |>
  filter(input_ng >= 49 & coverage.x >= 1 &
    path_block_manual_check == "pathology blocks match")

high_quality_results |>
  group_by(hrd_status_check) |>
  summarise(total = n())

# Results with quality control added
ggplot(
  high_quality_results,
  aes(x = myriad_gi_score, y = seqone_hrd_score)
) +
  geom_point(size = 3, alpha = 0.6, aes(colour = hrd_status_check)) +
  scale_colour_manual(values = c(safe_blue, safe_red)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = c(0, 25, 42, 50, 75, 100)
  ) +
  ylim(0, 1) +
  labs(
    x = "Myriad Genome Instability Score",
    y = "SeqOne HRD Score",
    title = "Comparison of Myriad vs SeqOne HRD Testing",
    subtitle = paste0(
      "Data for ", nrow(high_quality_results),
      " DNA inputs: path blocks match, >=49ng input, >=1x coverage"
    )
  ) +
  ggpubr::stat_cor(method = "pearson", label.x = 50, label.y = 0.1) +
  theme_bw() +
  theme(
    panel.grid = element_blank(), legend.position = "bottom",
    legend.title = element_blank()
  ) +
  guides(shape = guide_legend(ncol = 1)) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 42, linetype = "dashed")

# Coverage vs DNA input
ggplot(
  compare_results |>
    filter(path_block_manual_check == "pathology blocks match"),
  aes(x = input_ng, y = coverage.x)
) +
  geom_point(size = 3, alpha = 0.6, aes(
    shape = hrd_status_check,
    colour = hrd_status_check
  )) +
  scale_colour_manual(values = c("#CCCCCC", "#FF0000")) +
  theme_bw() +
  theme(
    panel.grid = element_blank(), legend.position = "bottom",
    legend.title = element_blank()
  ) +
  guides(shape = guide_legend(ncol = 1)) +
  ylim(0, 2)

# Total yield vs DNA input
ggplot(
  compare_results |>
    filter(path_block_manual_check == "pathology blocks match"),
  aes(x = input_ng, y = total_yield)
) +
  geom_point(size = 3, alpha = 0.6, aes(
    shape = hrd_status_check,
    colour = hrd_status_check
  )) +
  scale_colour_manual(values = c("#CCCCCC", "#FF0000")) +
  theme_bw() +
  theme(
    panel.grid = element_blank(), legend.position = "bottom",
    legend.title = element_blank()
  )

# DNA samples with lower than 50ng input
ggplot(
  compare_results |>
    filter(input_ng < 49 & !is.na(myriad_gi_score)),
  aes(x = myriad_gi_score, y = seqone_hrd_score)
) +
  geom_point(size = 3, alpha = 0.6, aes(
    shape = hrd_status_check,
    colour = hrd_status_check
  )) +
  scale_colour_manual(values = c("#CCCCCC", "#FF0000")) +
  theme_bw() +
  theme(
    panel.grid = element_blank(), legend.position = "bottom",
    legend.title = element_blank()
  ) +
  xlim(0, 100) +
  ylim(0, 1) +
  ggpubr::stat_cor(method = "pearson", label.x = 50, label.y = 0.25) +
  labs(title = "DNA samples with lower than 50ng input")

## Individual discrepant samples ----------------------------------------------------

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

coverage_plot <- compare_results |>
  filter(path_block_manual_check != "NA") |>
  ggplot(aes(x = reorder(shallow_sample_id, coverage.x), y = coverage.x)) +
  geom_point(size = 3, aes(
    shape = hrd_status_check,
    colour = hrd_status_check
  )) +
  scale_colour_manual(values = c("#CCCCCC", "#FF0000")) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(
    x = "Sample", y = "Coverage",
    caption = "Plot name: coverage_plot"
  ) +
  facet_wrap(~path_block_manual_check)

## Export tables --------------------------------------------------------------------

inconsistent_summary <- compare_results |>
  filter(hrd_status_check == "Seqone HRD status NOT consistent with Myriad") |>
  select(
    sample_id, worksheet, seqone_hrd_score, myriad_gi_score,
    seqone_hrd_status, myriad_hrd_status,
    path_block_manual_check, lga, lpc, ccne1, rad51b, coverage.x, percent_mapping,
    myriad_r_number, input_ng
  ) |>
  arrange(path_block_manual_check, sample_id)

export_timestamp(hrd_output_path, inconsistent_summary)

# Entire table of results
export_timestamp(hrd_output_path, compare_results)

# Remove patient identifiable information

results_patient_ids_removed <- compare_results |>
  select(-c(
    nhs_number, firstname, surname, myriad_patient_name,
    myriad_dob, myriad_filename
  ))

export_timestamp(hrd_output_path, results_patient_ids_removed)

## Impact of read length ------------------------------------------------------------

compare_results |>
  filter(downsampled == "No" & path_block_manual_check != "NA") |>
  ggplot(aes(x = reorder(shallow_sample_id, read_length), y = read_length)) +
  geom_point(size = 3, alpha = 0.5, aes(shape = hrd_status_check)) +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  labs(
    title = "SeqOne Read Length",
    x = "",
    y = "Read length"
  ) +
  facet_wrap(~path_block_manual_check)

## Comparison with SeqOne simplified model ------------------------------------------

simplified_model <- compare_results |>
  filter(!is.na(myriad_brca_status)) |>
  mutate(
    brca_status = case_when(
      myriad_brca_status == "POSITIVE" ~ "POSITIVE",
      myriad_brca_status == "NEGATIVE" ~ "NEGATIVE",
      # Seraseq controls and Biobank controls are BRCA negative
      myriad_brca_status == "" ~ "NEGATIVE"
    ),
    model_approximation = case_when(
      lga >= 18 ~ "POSITIVE",
      lga <= 14 ~ "NEGATIVE",
      lga < 18 & lga > 14 & lpc > 10 ~ "POSITIVE",
      lga < 18 & lga > 14 & lpc <= 10 ~ "NEGATIVE"
    ),
    approximation_consistent = ifelse(model_approximation == seqone_hrd_status,
      "Yes",
      "No"
    ),
    approximation_myriad_consistent = ifelse(model_approximation == myriad_hrd_status,
      "Yes",
      "No"
    )
  ) |>
  filter(brca_status == "NEGATIVE") |>
  mutate(ccne1_rounded = round(ccne1, 0))

simplified_model_summary <- simplified_model |>
  group_by(approximation_consistent) |>
  summarise(total = n())

lga_vs_lpc_plot <- ggplot(simplified_model, aes(
  x = lga,
  y = lpc
)) +
  geom_point(aes(colour = seqone_hrd_status, size = ccne1_rounded),
    alpha = 0.6
  ) +
  scale_colour_manual(values = c(safe_blue, safe_red)) +
  ylim(0, 40) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(
    x = "Large Genomic Alterations",
    y = "Loss of Parental Copy",
    # title = "Comparison of Seqone model approximation to pipeline output",
    # subtitle = paste0("Data for ", nrow(simplified_model), " BRCA negative DNA inputs"),
    colour = "SeqOne pipeline status",
    size = "CCNE1"
  ) +
  geom_segment(aes(x = 18, y = 0, xend = 18, yend = 10), linetype = "dashed") +
  geom_segment(aes(x = 15, y = 10, xend = 15, yend = 40), linetype = "dashed") +
  geom_segment(aes(x = 15, y = 10, xend = 18, yend = 10), linetype = "dashed") +
  scale_size_binned(breaks = c(1, 2, 3, 4, 5, 10))

save_hrd_plot(lga_vs_lpc_plot, input_width = 15.5, input_height = 15)

# Does the approximation agree with Myriad HRD status?

myriad_v_approximation_plot <- ggplot(simplified_model |>
  filter(path_block_manual_check == "pathology blocks match" &
    !is.na(myriad_hrd_status)), aes(x = lga, y = lpc)) +
  geom_point(aes(colour = myriad_hrd_status),
    alpha = 0.6, size = 3
  ) +
  scale_colour_manual(values = c(safe_blue, safe_red)) +
  ylim(0, 40) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(
    x = "Large Genomic Alterations",
    y = "Loss of Parental Copy",
    colour = "Myriad HRD status",
    title = "SeqOne model approximation versus Myriad HRD status"
  ) +
  geom_segment(aes(x = 18, y = 0, xend = 18, yend = 10), linetype = "dashed") +
  geom_segment(aes(x = 15, y = 10, xend = 15, yend = 40), linetype = "dashed") +
  geom_segment(aes(x = 15, y = 10, xend = 18, yend = 10), linetype = "dashed")

# Does the approximation agree with SeqOne HRD status?

seqon_v_approximation_plot <- ggplot(simplified_model |>
  filter(path_block_manual_check == "pathology blocks match" &
    !is.na(myriad_hrd_status)), aes(x = lga, y = lpc)) +
  geom_point(aes(colour = seqone_hrd_status),
    alpha = 0.6, size = 3
  ) +
  scale_colour_manual(values = c(safe_blue, safe_red)) +
  ylim(0, 40) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(
    x = "Large Genomic Alterations",
    y = "Loss of Parental Copy",
    colour = "SeqOne HRD status",
    title = "SeqOne model approximation versus SeqOne HRD status"
  ) +
  geom_segment(aes(x = 18, y = 0, xend = 18, yend = 10), linetype = "dashed") +
  geom_segment(aes(x = 15, y = 10, xend = 15, yend = 40), linetype = "dashed") +
  geom_segment(aes(x = 15, y = 10, xend = 18, yend = 10), linetype = "dashed")

merge_plot <- ggpubr::ggarrange(
  seqon_v_approximation_plot,
  myriad_v_approximation_plot
)

## Seraseq controls -----------------------------------------------------------------

seraseq_control_data <- compare_results |>
  filter(control_type == "Seraseq control") |>
  mutate(firstname_factor = factor(firstname, levels = c(
    "FFPE HRD Negative",
    "Low-Positive FFPE HRD",
    "High-Positive FFPE HRD"
  )))

ggplot(seraseq_control_data, aes(x = seqone_hrd_score, y = insert_size)) +
  geom_point(size = 3, aes(
    shape = hrd_status_check,
    colour = worksheet
  )) +
  facet_wrap(~firstname_factor) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(title = "Seraseq Controls: repeat SeqOne data") +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  xlim(0, 1)

## Biobank controls -----------------------------------------------------------------

compare_results |>
  filter(surname %in% biobank_names) |>
  ggplot(aes(x = surname, y = seqone_hrd_score)) +
  geom_point(size = 3) +
  ylim(0, 1) +
  labs(x = "", y = "SeqOne HRD score") +
  theme_bw()

## Myriad results -------------------------------------------------------------------

collated_myriad_info_mod |>
  ggplot(aes(
    x = reorder(nhs_number, myriad_gi_score),
    y = myriad_gi_score
  )) +
  geom_point(size = 3, aes(shape = myriad_hrd_status)) +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  geom_hline(yintercept = 42, linetype = "dashed")

## QC metrics -----------------------------------------------------------------------

# Million reads
ggplot(compare_results, aes(
  x = reorder(shallow_sample_id, million_reads),
  y = million_reads
)) +
  geom_point(size = 3, aes(colour = hrd_status_check)) +
  labs(x = "", y = "million_reads") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90),
    legend.position = "bottom",
    legend.title = element_blank()
  )

# Read length
ggplot(compare_results, aes(
  x = reorder(shallow_sample_id, read_length),
  y = read_length
)) +
  geom_point(size = 3, aes(colour = hrd_status_check)) +
  # scale_colour_manual(values = c(negative_colour,
  # positive_colour)) +
  labs(x = "", y = "read_length") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90),
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  facet_wrap(~worksheet)

# Insert size
ggplot(compare_results, aes(
  x = reorder(shallow_sample_id, insert_size),
  y = insert_size
)) +
  geom_point(size = 3, aes(colour = hrd_status_check)) +
  labs(x = "", y = "insert_size") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90),
    legend.position = "bottom",
    legend.title = element_blank()
  )

# Percent Q30
ggplot(compare_results, aes(
  x = reorder(shallow_sample_id, percent_q30),
  y = percent_q30
)) +
  geom_point(size = 3, aes(colour = hrd_status_check)) +
  labs(x = "", y = "percent_q30") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90),
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  ylim(0, 100)

# Percent aligned
ggplot(compare_results, aes(
  x = reorder(shallow_sample_id, percent_aligned),
  y = percent_aligned
)) +
  geom_point(size = 3, aes(colour = hrd_status_check)) +
  labs(x = "", y = "percent_aligned") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90),
    legend.position = "bottom",
    legend.title = element_blank()
  )

# Percent dups
ggplot(compare_results, aes(
  x = reorder(shallow_sample_id, percent_dups),
  y = percent_dups
)) +
  geom_point(size = 3, aes(colour = hrd_status_check)) +
  labs(x = "", y = "percent_dups") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90),
    legend.position = "bottom",
    legend.title = element_blank()
  )

# Coverage
ggplot(compare_results, aes(
  x = reorder(shallow_sample_id, coverage.y),
  y = coverage.y
)) +
  geom_point(size = 3, aes(colour = hrd_status_check)) +
  labs(x = "", y = "coverage.y") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90),
    legend.position = "bottom",
    legend.title = element_blank()
  )

ggplot(compare_results, aes(
  x = coverage.x,
  y = read_length
)) +
  geom_point(size = 3, aes(colour = hrd_status_check)) +
  labs(x = "coverage", y = "read_length") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90),
    legend.position = "bottom",
    legend.title = element_blank()
  )

## Sample extraction information ----------------------------------------------------

# Get tissues
seqone_tissue_types <- unique(seqone_dlms_info$tissue)

tissue_query <- sqlQuery(
  channel = moldb_connection,
  query = paste0(
    "SELECT * FROM MolecularDB.dbo.TissueTypes WHERE TissueTypeId IN (",
    paste(seqone_tissue_types, collapse = ", "),
    ")"
  )
) |>
  janitor::clean_names()

# Get extraction batch IDs
extraction_batch_id_table <- sqlQuery(
  channel = moldb_connection,
  query = paste0(
    "SELECT * FROM MolecularDB.dbo.MOL_Extractions WHERE LABNO IN (",
    paste(seqone_dlms_info$dlms_dna_number, collapse = ", "),
    ")"
  )
) |>
  arrange(LabNo) |>
  janitor::clean_names() |>
  dplyr::rename(
    dlms_dna_number = lab_no,
    extraction_batch_id = extraction_batch_fk
  )


extraction_batches <- unique(extraction_batch_id_table$extraction_batch_id)

# Get extraction batch dates
extraction_batch_date_table <- sqlQuery(
  channel = moldb_connection,
  query = paste0(
    "SELECT * FROM MolecularDB.dbo.MOL_ExtractionBatches WHERE ExtractionBatchId IN (",
    paste(extraction_batches, collapse = ", "),
    ")"
  )
) |>
  janitor::clean_names()

extraction_batch_info <- extraction_batch_id_table |>
  left_join(extraction_batch_date_table, by = "extraction_batch_id") |>
  # Get only extraction batch IDs used for Cobas extractions
  filter(extraction_method_fk == 25)

# 1 DNA number (23033279) is on 2 extraction batches

seqone_dlms_extractions <- seqone_dlms_info |>
  left_join(tissue_query |>
    dplyr::rename(tissue = tissue_type_id), by = "tissue") |>
  left_join(extraction_batch_info, by = "dlms_dna_number")

# Plot of tissue types
ggplot(seqone_dlms_extractions, aes(x = tissue_type, y = )) +
  geom_bar()

ggplot(seqone_dlms_extractions, aes(x = run_date, y = concentration)) +
  geom_point()

## Impact of DNA concentration ------------------------------------------------------

ggplot(
  compare_results |>
    filter(path_block_manual_check == "pathology blocks match"),
  aes(x = input_ng, y = coverage.x)
) +
  geom_point(size = 3, aes(
    shape = hrd_status_check,
    colour = hrd_status_check
  )) +
  scale_colour_manual(values = c("#CCCCCC", "#FF0000")) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  ylim(0, 2) +
  xlim(0, 55) +
  facet_wrap(~path_block_manual_check)
# geom_point(data=compare_results[compare_results$dlms_dna_number == 23034142,],
# aes(input_ng, coverage.x),
# =21, fill=NA, size=5, colour="red", stroke=2)

ggplot(
  compare_results |>
    filter(path_block_manual_check == "pathology blocks match"),
  aes(x = input_genomes, y = lga)
) +
  geom_point(size = 3, aes(
    shape = hrd_status_check,
    colour = hrd_status_check
  )) +
  scale_colour_manual(values = c("#CCCCCC", "#FF0000")) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  facet_wrap(~path_block_manual_check)

ggplot(
  compare_results |>
    filter(path_block_manual_check == "pathology blocks match"),
  aes(
    x = reorder(dlms_dna_number, qubit_dna_ul),
    y = qubit_dna_ul
  )
) +
  geom_point(size = 3, aes(
    shape = hrd_status_check,
    colour = hrd_status_check
  )) +
  scale_colour_manual(values = c("#CCCCCC", "#FF0000")) +
  geom_hline(yintercept = 3.3, linetype = "dashed") +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  labs(x = "")

## Tumour BRCA referral DNA concentrations ------------------------------------------

brca_query <- "SELECT * FROM MolecularDB.dbo.Samples WHERE DISEASE IN (204)"

tbrca_data <- sqlQuery(
  channel = moldb_connection,
  query = brca_query
) |>
  janitor::clean_names()

dna_qc_threshold <- round(50 / 15, 0)

tbrca_data_mod <- tbrca_data |>
  filter(!is.na(concentration)) |>
  mutate(pass_qc = ifelse(concentration >= dna_qc_threshold, "Yes", "No"))

samples_passing_qc <- nrow(tbrca_data_mod[tbrca_data_mod$pass_qc == "Yes", ])

samples_failing_qc <- nrow(tbrca_data_mod[tbrca_data_mod$pass_qc == "No", ])

fail_rate <- round((samples_failing_qc /
  (samples_passing_qc + samples_failing_qc)) * 100, 1)

ggplot(tbrca_data_mod, aes(x = disease, y = concentration)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  labs(
    y = "DNA concentration (ng/ul)",
    x = "Tumour BRCA referrals",
    title = paste0(
      "DNA concentrations for ",
      nrow(tbrca_data_mod), " tumour BRCA referrals"
    ),
    subtitle = paste0(
      "Samples below ", dna_qc_threshold, " ng/ul: ", samples_failing_qc,
      " (", fail_rate, "%)"
    ),
    caption = paste0(
      "Median DNA concentration: ", median(tbrca_data_mod$concentration),
      " ng/ul"
    )
  ) +
  ylim(0, 700) +
  geom_hline(yintercept = 3.3, linetype = "dashed")

min(tbrca_data_mod$date_in)

## Spread of Myriad GI scores -------------------------------------------------------

tbrca_data_collection <- read_excel(
  paste0(
    hrd_data_path,
    "HRD TBRCA data collection Manchester_NEW_from Oct2022_2023.xlsx"
  ),
  skip = 1
) |>
  janitor::clean_names()

tbrca_data_collection_clean <- tbrca_data_collection |>
  filter(!gis_score_numerical_value_or_fail_not_tested %in% c(
    "Fail", "Inconclusive",
    "Not tested", NA
  )) |>
  filter(t_brca_mutation_status != "Fail") |>
  mutate(
    gi_score = as.numeric(gis_score_numerical_value_or_fail_not_tested),
    brca_mutation_clean = case_when(
      t_brca_mutation_status %in%
        c("Pathogenic BRCA1", "Pathogenic BRCA2") ~ "BRCA positive",
      t_brca_mutation_status %in%
        c("No mutation detected", "no mutation detected") ~ "BRCA negative"
    ),
    borderline = ifelse(gi_score >= 37 & gi_score <= 47, "Yes", "No")
  )

tbrca_data_collection_clean |>
  group_by(borderline) |>
  summarise(total = n())

243 / (sum(1864, 243))

myriad_gi_profile_plot <- tbrca_data_collection_clean |>
  ggplot(aes(x = gi_score, y = )) +
  geom_histogram(binwidth = 1, aes(fill = gis_pos_neg)) +
  scale_fill_manual(values = c(safe_blue, safe_red)) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = c(0, 25, 42, 50, 75, 100)
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom"
  ) +
  # geom_vline(xintercept = 42, linetype = "dashed") +
  labs(
    y = "Number of samples",
    x = "Myriad GI score",
    caption = paste0("Data for ", nrow(tbrca_data_collection_clean), " samples shown"),
    fill = "Myriad GI Status",
    title = "Myriad GI Scores for North West GLH Samples"
  )

save_hrd_plot(myriad_gi_profile_plot)

compare_results |>
  filter(!is.na(seqone_hrd_score)) |>
  group_by(seqone_hrd_score, seqone_hrd_status) |>
  summarise(total = n()) |>
  ggplot(aes(x = seqone_hrd_score, y = total)) +
  geom_col(aes(fill = seqone_hrd_status), width = 0.01) +
  # geom_smooth() +
  scale_fill_manual(values = c(safe_blue, safe_red)) +
  scale_x_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom"
  ) +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  labs(
    y = "Number of samples",
    x = "Seqone HRD score",
    fill = "SeqOne GI Status",
    title = "Seqone HRD Scores for North West GLH Samples"
  )

compare_results |>
  filter(!is.na(seqone_hrd_score)) |>
  mutate(borderline = ifelse(seqone_hrd_score >= 0.4 & seqone_hrd_score <= 0.6,
    "Yes", "No"
  )) |>
  group_by(borderline) |>
  summarise(total = n())

11 / (11 + 78)

## Coverage over worksheets ---------------------------------------------------------

ggplot(compare_results, aes(x = worksheet, y = coverage.x)) +
  geom_boxplot() +
  ylim(0, 2)

## RAD51B and CCNE1 -----------------------------------------------------------------

ggplot(collated_seqone_info, aes(x = rad51b, y = ccne1)) +
  geom_point() +
  theme_bw() +
  xlim(0, 10) +
  ylim(0, 10)

ggplot(collated_seqone_info |>
  filter(ccne1 < 5), aes(
  x = ,
  y = ccne1
)) +
  geom_boxplot() +
  scale_y_continuous(breaks = c(0:5))

ggplot(collated_seqone_info, aes(
  x = ccne1,
  y = seqone_hrd_score
)) +
  geom_point(size = 3, alpha = 0.6, aes(colour = seqone_hrd_status)) +
  scale_colour_manual(values = c(safe_blue, safe_red)) +
  theme_bw() +
  scale_x_continuous(breaks = c(0:10))

## SeqOne amended results -----------------------------------------------------------

seqone_myriad_plot <- compare_results |>
  filter(path_block_manual_check == "pathology blocks match") |>
  ggplot(aes(seqone_hrd_status, myriad_hrd_status)) +
  geom_jitter(
    size = 3, alpha = 0.6,
    width = 0.2, height = 0.2
  ) +
  theme_bw() +
  labs(
    title = "Comparison of SeqOne original results with Myriad",
    subtitle = "9 discrepant results"
  )

seqone_amended_myriad_plot <- compare_results |>
  filter(path_block_manual_check == "pathology blocks match") |>
  ggplot(aes(seqone_hrd_status_amended, myriad_hrd_status)) +
  geom_jitter(
    size = 3, alpha = 0.6,
    width = 0.2, height = 0.2
  ) +
  theme_bw() +
  labs(
    title = "Comparison of SeqOne amended results with Myriad",
    subtitle = "False positive is sample 23016526: 2.2ng/ul and 0.47x coverage"
  )

ggarrange(seqone_myriad_plot, seqone_amended_myriad_plot, nrow = 1)

line_df <- data.frame(
  x = c(18, 17, 16, 15, 14, 13, 18, 17, 16, 15, 14),
  y = c(0, 4, 9, 13, 18, 23, 4, 9, 13, 18, 23),
  xend = c(18, 17, 16, 15, 14, 13, 17, 16, 15, 14, 13),
  yend = c(4, 9, 13, 18, 23, 28, 4, 9, 13, 18, 23)
)

lga_vs_lpc <- ggplot(compare_results, aes(lga, lpc)) +
  geom_point(aes(colour = seqone_hrd_status),
    size = 3
  ) +
  scale_colour_manual(values = c(safe_blue, safe_red)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  geom_segment(
    data = line_df,
    mapping = aes(x = x, y = y, xend = xend, yend = yend),
    linetype = "dashed"
  ) +
  labs(
    title = "Original pipeline",
    subtitle = "CCNE1 and RAD51B included"
  )


lga_vs_lpc_amended <- ggplot(compare_results, aes(lga_amended, lpc_amended)) +
  geom_point(aes(colour = seqone_hrd_status_amended),
    size = 3
  ) +
  scale_colour_manual(values = c(safe_blue, "#CCCCCC", safe_red)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  geom_segment(
    data = line_df,
    mapping = aes(x = x, y = y, xend = xend, yend = yend),
    linetype = "dashed"
  ) +
  labs(
    title = "New pipeline",
    subtitle = "CCNE1 and RAD51B removed, QC step added"
  )

ggarrange(lga_vs_lpc, lga_vs_lpc_amended,
  nrow = 1
)

lga_plot <- ggplot(compare_results, aes(lga_amended, lga)) +
  geom_point(size = 3, aes(colour = seqone_hrd_status_amended)) +
  scale_colour_manual(values = c(safe_blue, "#CCCCCC", safe_red)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(title = "Large genomic alterations") +
  ylim(0, 45) +
  xlim(0, 45)

lpc_plot <- ggplot(compare_results, aes(lpc_amended, lpc)) +
  geom_point(size = 3, aes(colour = seqone_hrd_status_amended)) +
  scale_colour_manual(values = c(safe_blue, "#CCCCCC", safe_red)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(title = "Loss of parental copy") +
  ylim(0, 45) +
  xlim(0, 45)

ggarrange(lga_plot, lpc_plot, nrow = 1)

compare_results |>
  filter(seqone_hrd_status_amended == "NON-CONCLUSIVE") |>
  select(
    shallow_sample_id, coverage.x, input_ng, seqone_hrd_status_amended,
    low_tumor_fraction_returns_a_warning_only
  ) |>
  view()

compare_results |>
  filter(low_tumor_fraction_returns_a_warning_only == "YES") |>
  select(
    shallow_sample_id, coverage.x, input_ng, seqone_hrd_status_amended,
    low_tumor_fraction_returns_a_warning_only
  ) |>
  view()

## Sensitivity and specificity ------------------------------------------------------

perform_sensitivity_calcs(compare_results |>
  filter(path_block_manual_check == "pathology blocks match"))

perform_sensitivity_calcs(compare_results |>
  filter(path_block_manual_check == "pathology blocks match" &
    input_ng > 49 &
    coverage.x >= 0.5))
