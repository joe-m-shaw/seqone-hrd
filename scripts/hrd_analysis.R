# Homologous Recombination Deficiency Validation: Analysis
# joseph.shaw3@nhs.net

# Setup -----------------------------------------------------------------------------

rm(list = ls())

## Functions ------------------------------------------------------------------------

source("functions/hrd_functions.R")

## Packages and filepaths -----------------------------------------------------------

library("ggpubr")
library("readxl")
library("RODBC")

## DLMS connection ------------------------------------------------------------------

# Connection setup for PC38698
moldb_connection <- RODBC::odbcConnect(dsn = "moldb")

# Collate data ----------------------------------------------------------------------

## Collate Myriad data --------------------------------------------------------------

myriad_report_files <- list.files(myriad_reports_location, full.names = TRUE)

collated_myriad_info <- myriad_report_files |>
  map(read_myriad_report) |>
  list_rbind()

# Check all files included
stopifnot(setdiff(
  basename(myriad_report_files),
  collated_myriad_info$myriad_filename
) == 0)

## Collate SeqOne data --------------------------------------------------------------

seqone_report_files <- list.files(seqone_report_location, full.names = TRUE)

collated_seqone_info <- seqone_report_files |>
  map(read_seqone_report) |>
  list_rbind()

# Check all files collated
stopifnot(setdiff(
  basename(seqone_report_files),
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
  collated_seqone_info$dlms_dna_number
) == 0)

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

# export_timestamp(hrd_data_path, export_for_check)

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
  extract_kapa_data("135001", 31),
  extract_kapa_data("135498", 7)
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

join_tables <- seqone_mod |>
  left_join(collated_myriad_info_mod, by = "nhs_number") |>
  left_join(path_block_check, by = "dlms_dna_number") |>
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
  ) 

# Analyse Data ----------------------------------------------------------------------

## Compare SeqOne and Myriad results ------------------------------------------------

compare_results <- join_tables |> 
  filter(downsampled == "No") |>
  mutate(
    
    # Outome of results from original pipeline
    outcome_v1 = case_when(
      
      seqone_hrd_status == pos_text & myriad_hrd_status == pos_text ~true_pos_text,
      
      seqone_hrd_status == neg_text & myriad_hrd_status == neg_text ~true_neg_text,
      
      seqone_hrd_status == pos_text & myriad_hrd_status == neg_text ~false_pos_text,
      
      seqone_hrd_status == neg_text & myriad_hrd_status == pos_text ~false_neg_text,
      
      TRUE ~ "other"),
    
    outcome_v2 = case_when(
      
      seqone_hrd_status_amended == pos_text & myriad_hrd_status == pos_text ~true_pos_text,
      
      seqone_hrd_status_amended == neg_text & myriad_hrd_status == neg_text ~true_neg_text,
      
      seqone_hrd_status_amended == pos_text & myriad_hrd_status == neg_text ~false_pos_text,
      
      seqone_hrd_status_amended == neg_text & myriad_hrd_status == pos_text ~false_neg_text,
      
      seqone_hrd_status_amended == incon_text & myriad_hrd_status == pos_text ~incon_pos_text,
      
      seqone_hrd_status_amended == incon_text & myriad_hrd_status == neg_text ~incon_neg_text,
      
      TRUE ~ "other")
    
    )
      

## Sensitivity and specificity ------------------------------------------------------

coverage_threshold <- 0.5

input_threshold <- 47

v1_results <- compare_results |> 
  filter(path_block_manual_check == "pathology blocks match") |> 
  compare_tests(outcome_v1)

v1_results_filtered <- compare_results |> 
  filter(path_block_manual_check == "pathology blocks match" &
           coverage.x >= coverage_threshold & input_ng >= input_threshold) |> 
  compare_tests(outcome_v1)

v2_results <- compare_results |> 
  filter(path_block_manual_check == "pathology blocks match") |> 
  compare_tests(outcome_v2)

v2_results_filtered <- compare_results |> 
  filter(path_block_manual_check == "pathology blocks match" &
           coverage.x >= coverage_threshold & input_ng >= input_threshold) |> 
  compare_tests(outcome_v2)

add_version <- function(input_table, version_text) {
  
  input_table |> 
    mutate(Analysis = version_text) |> 
    relocate(Analysis)
  
}

metric_table <- rbind(add_version(v1_results[[2]], "SomaHRD v1"),
                      add_version(v1_results_filtered[[2]], "SomaHRD v1, thresholds applied"),
                      add_version(v2_results[[2]], "SomaHRD v2"),
                      add_version(v2_results_filtered[[2]], "SomaHRD v2, thresholds applied"))

# export_timestamp(hrd_output_path, metric_table)

# export_timestamp(hrd_output_path, v1_results[[1]])

# export_timestamp(hrd_output_path, v2_results[[1]])

## Myriad and SeqOne score correlation ----------------------------------------------

compare_results |> 
  filter(path_block_manual_check == "pathology blocks match" & 
           coverage.x >= coverage_threshold & input_ng >= input_threshold) |> 
  ggplot(aes(myriad_gi_score, seqone_hrd_score)) +
  geom_point(size = 3, alpha = 0.7) +
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
    title = "Comparison of Myriad vs SeqOne HRD Scores",
    subtitle = str_c("DNA inputs with coverage >= ", coverage_threshold,
                     "X; DNA input >= ", input_threshold, "ng")) +
  ggpubr::stat_cor(method = "pearson", label.x = 75, label.y = 0.1)

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
    alpha = 0.5, size = 2,
    aes(
      shape = seqone_hrd_status
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

#save_hrd_plot(repeat_facet_plot, input_width = 15, input_height = 15)

## Intra-run variation --------------------------------------------------------------

intra_run_table <- compare_results |> 
  filter(worksheet == "WS133557") |> 
  filter(base::duplicated(dlms_dna_number, fromLast = TRUE) |
           base::duplicated(dlms_dna_number, fromLast = FALSE)) |> 
  select(sample_id, seqone_hrd_status, seqone_hrd_score, lga, lpc, coverage.x) |> 
  arrange(sample_id)

export_timestamp(input = intra_run_table)

## Comparison with SeqOne simplified model ------------------------------------------

## Seraseq controls -----------------------------------------------------------------

seraseq_control_data <- compare_results |>
  filter(control_type == "Seraseq control") |>
  mutate(firstname_factor = factor(firstname, levels = c(
    "FFPE HRD Negative",
    "Low-Positive FFPE HRD",
    "High-Positive FFPE HRD"
  )))

ggplot(seraseq_control_data, aes(x = worksheet, y = seqone_hrd_score)) +
  geom_point(size = 3) +
  facet_wrap(~firstname_factor) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(title = "Seraseq Controls: repeat SeqOne data")

## Biobank controls -----------------------------------------------------------------

compare_results |>
  filter(surname %in% biobank_names) |>
  ggplot(aes(x = surname, y = seqone_hrd_score)) +
  geom_point(size = 3) +
  ylim(0, 1) +
  labs(x = "", y = "SeqOne HRD score") +
  theme_bw()

## QC metrics -----------------------------------------------------------------------

# Add function to plot results by QC metric

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

#save_hrd_plot(myriad_gi_profile_plot)

## SeqOne amended results -----------------------------------------------------------

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
