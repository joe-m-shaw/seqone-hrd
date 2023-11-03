# HRD Project Functions

library(pdftools)
library(tidyverse)
library(assertthat)

# Standard text ---------------------------------------------------------------------

pos_text <- "POSITIVE"

neg_text <- "NEGATIVE"

incon_text <- "NON-CONCLUSIVE"

true_pos_text <- "true_positive"

true_neg_text <- "true_negative"

false_pos_text <- "false_positive"

false_neg_text <- "false_negative"

incon_pos_text <- "inconclusive_positive"

incon_neg_text <- "inconclusive_negative"

consistent_text <- "Seqone HRD status consistent with Myriad"

inconsistent_text <- "Seqone HRD status NOT consistent with Myriad"

inconclusive_text <- "SeqOne HRD status inconclusive"


# Filepaths -------------------------------------------------------------------------

dev_team_path <- "S:/central shared/Genetics/Mol_Shared/Development.Team/"

seqone_folder <- paste0(dev_team_path, 
                        "SeqOne Homologous Recombination Deficiency Validation/")

hrd_project_path <- paste0(seqone_folder, "HRD R script files/")

hrd_data_path <- paste0(hrd_project_path, "data/")

hrd_output_path <- paste0(hrd_project_path, "outputs/")

hrd_plot_path <- paste0(hrd_project_path, "plots/")

myriad_reports_location <- paste0(hrd_data_path, "myriad_reports/")

seqone_report_location <- paste0(hrd_data_path, "seqone_reports/")

# CSV timestamp ---------------------------------------------------------------------

export_timestamp <- function(filepath = hrd_output_path, input) {
  write.csv(input,
    file = paste0(
      filepath,
      format(Sys.time(), "%Y_%m_%d_%H_%M_%S"),
      "_",
      deparse(substitute(input)), ".csv"
    ),
    row.names = FALSE
  )
}


# Reading Myriad PDFs ---------------------------------------------------------------

check_na <- function(input_table) {
  assert_that(sum(is.na(input_table)) == 0,
    msg = paste0("NA values present: file ", file)
  )
}

read_myriad_report <- function(file) {
  # Use pdf_text to read PDF as a single string per page

  myriad_report_text <- pdftools::pdf_text(pdf = file)

  page1 <- myriad_report_text[[1]]

  page2 <- myriad_report_text[[2]]

  page3 <- myriad_report_text[[3]]

  # iGene R Number
  
  r_number_regex <- regex(
    r"[
    Patient\sID:        # Escape space with backslash s
    \s{6}               # 6 spaces
    (R\d{2}-\w{4})      # R number grouped
    ]",
    comments = TRUE
  )
  
  myriad_r_number <- str_extract(page2, r_number_regex, group = 1)

  # NHS Number
  
  nhs_no_regex <- regex(
    r"[
    NHS\sNo:        
    .{10,20}                    # Size of whitespace between No: and number can vary
                                # Use . instead of backslash s as samples can have a leading
                                # 0 before NHS number. Example: R22-031L
    (\d{3}\s\d{3}\s\d{4})       # Grouped NHS number format
    ]",
    comments = TRUE
  )
  
  nhs_no_char <- str_extract(page1, nhs_no_regex, group = 1)
  
  nhs_no_double <- parse_number(nhs_no_char, locale = locale(grouping_mark = " "))
  
  # Patient name
  
  patient_name_regex <- regex(
    r"[
    Patient\sName:
    \s{10,15}           # Variable whitespace
    (\D{5,25})          # Names as non-digits with variable lengths
    \n                  # New line marks end of the name
    ]",
    comments = TRUE
  )
  
  myriad_patient_name <- str_extract(page1, patient_name_regex, group = 1)

  # Patient date of birth
  
  dob_regex <- regex(
    r"[
    Date\sof\sBirth:
    \s{10,15}                 # Variable whitespace
    (\d{2}/\d{2}/\d{4})       # DOB in dd/mm/yyyy format
    \n                        # Mark end of DOB with newline
    ]",
    comments = TRUE
  )
  
  myriad_dob <- str_extract(page1, dob_regex, group = 1)

  # Pathology Block
  
  path_block_1_regex <- regex(
    r"[
    Pathology\sNo:
    \s{10,30}         # Variable whitespace
    (.{5,25})         # Path blocks can have various formats from different labs
    \n                # Specify end of block number with new line
    ]",
    comments = TRUE
  )
  
  myriad_pathology_block_pg1 <- str_extract(page1, path_block_1_regex, group = 1)
  
  path_block_2_regex <- regex(
    r"[
    (Block\(s\)|Specimen\(s\))  # Some say Block(s) whilst others say Specimen(s)
    \sAnalyzed:\s
    (.{5,20})                   # Path blocks can have various formats from different labs
    \n                          # Specify end of block number with new line
    ]",
    comments = TRUE
  )
  
  myriad_pathology_block_pg2 <- str_extract(page2, path_block_2_regex, group = 2)

  # GI score
  
  gi_score_regex <- regex(
    r"[
    Patient\sGenomic\sInstability\sScore:
    \s
    (\d{1,2})                               # Score can be 1 digit (i.e. 7) or 2 (73)
    \n                                      # Mark end of score with newline
    ]",
    comments = TRUE
  )
  
  gi_score_char <- str_extract(str_c(page2, page3), gi_score_regex, group = 1)
  
  gi_score_double <- parse_number(gi_score_char)

  assert_that(is.na(gi_score_double) == FALSE,
    msg = paste0("GI score is NA: file ", basename(file))
  )

  assert_that(gi_score_double >= 0, gi_score_double <= 100, msg = paste0(
    "GI score outside 0-100 range: file ", basename(file)
  ))

  # HRD status
  
  hrd_status_regex <- regex(
    r"[
    Myriad\sHRD\sStatus:\s
    (NEGATIVE | POSITIVE)
    ]",
    comments = TRUE
  )
  
  myriad_hrd_status <- str_extract(page2, hrd_status_regex, group = 1)

  # BRCA status
  
  brca_status_regex <- regex(
    r"[
    Tumor\sMutation\sBRCA1/BRCA2\sStatus:\s
    (NEGATIVE | POSITIVE)
    ]",
    comments = TRUE
  )
  
  myriad_brca_status <- str_extract(page2, brca_status_regex, group = 1)

  output <- data.frame(
    "myriad_r_number" = myriad_r_number,
    "myriad_patient_name" = myriad_patient_name,
    "myriad_dob" = myriad_dob,
    "nhs_number" = nhs_no_double,
    "myriad_pathology_block_pg1" = myriad_pathology_block_pg1,
    "myriad_pathology_block_pg2" = myriad_pathology_block_pg2,
    "myriad_gi_score" = gi_score_double,
    "myriad_hrd_status" = myriad_hrd_status,
    "myriad_brca_status" = myriad_brca_status,
    "myriad_filename" = basename(file)
  )

  check_na(output)

  return(output)
}


# Reading SeqOne PDFs ---------------------------------------------------------------

check_version <- function(version) {
  
  assert_that(version %in% c("1.1", "1.2"))
  
}

get_hrd_score <- function(page, version) {
  
  check_version(version)
  
  hrd_score_regex_1_1 <- regex(
    r"[
    CLASS\n
    \s{78,83}                     # Variable whitespace
    ((0\.\d{1,2})|(\d{1}))        # Variable formats: 0.99, 0.9, 1
                                  # Use backslash . to specify decimal point
    ]",
    comments = TRUE
  )

  hrd_score_regex_1_2 <- regex(
    r"[
    HRD\sSummary\n
    \s{67}
    (\d{1}\.\d{2})
    ]",
    comments = TRUE
  )
  
  if(version == "1.1") {
    
    hrd_score_regex <- hrd_score_regex_1_1
    
  }
  
  if(version == "1.2") {
    
    hrd_score_regex <- hrd_score_regex_1_2
    
  }
  
  hrd_score_char <- str_extract(page, hrd_score_regex, group = 1)
  
  hrd_score_double <- parse_number(hrd_score_char, locale = locale(decimal_mark = "."))

  return(hrd_score_double)
  
}

get_hrd_status <- function(page) {
  
  hrd_status_regex <- regex(
    r"[
    SeqOne\sHRD\sStatus1\s:\s
    (NEGATIVE | POSITIVE)
    ]",
    comments = TRUE
  )
  
  seqone_hrd_status <- str_extract(page, hrd_status_regex, group = 1)
  
  return(seqone_hrd_status)
  
}

get_lga <- function(page) {
  
  lga_regex <- regex(
    r"[
    LGA\sStatus
    \s{35, 38}      # Variable whitespace - 35 in v1.2
    (\d{1,2})       # Group LGA status
    ]",
    comments = TRUE
  )
  
  lga <- parse_number(str_extract(page, lga_regex, group = 1))
  
  return(lga)
  
}

get_lpc <- function(page) {
  
  lpc_regex <- regex(
    r"[
    LPC\sStatus
    \s{35,38}            # Variable whitespace between versions
    (\d{1,2})
    ]",
    comments = TRUE
  )
  
  lpc <- parse_number(str_extract(page, lpc_regex, group = 1))
  
  return(lpc)
  
}

get_ncc <- function(page) {
  
  ncc_regex <- regex(
    r"[
    %\sof\stumoral\scells
    \s{23,44}               # Variable whitespace between versions
    (\d{2})                 # Grouped NCC value
    %
    ]",
    comments = TRUE
  )
  
  seqone_ncc <- parse_number(str_extract(page, ncc_regex, group = 1))
  
  assert_that(seqone_ncc >= 20, seqone_ncc <= 100)
  
  return(seqone_ncc)
  
}

get_coverage <- function(page) {
  
  coverage_regex <- regex(
    r"[
    Coverage
    \s{33,54}                      # Variable whitespace between versions
    ((\d{1}\.\d{1,2}) | (\d{1}))   # Variable format 1.57, 1.5, 1
    X
    ]",
    comments = TRUE
  )
  
  coverage <- parse_number(str_extract(page, coverage_regex, group = 1),
                           locale = locale(decimal_mark = "."))
  
  return(coverage)
  
}

get_percent_mapping <- function(page) {
  
  percent_map_regex <- regex(
    r"[
    %\scorrect\smapping
    \s{24,45}                     # Variable whitespace between versions
    ((\d{2}\.\d{1})|(\d{2}))      # Format 97.2, 97
                                  # Assume percent mapping always above 10
    %
    ]",
    comments = TRUE
  )
  
  percent_mapping <- parse_number(str_extract(page, percent_map_regex, group = 1),
                                  locale = locale(decimal_mark = "."))
  
  assert_that(percent_mapping >= 0, percent_mapping <= 100)
  
  return(percent_mapping)
  
}

get_shallow_sample_id <- function(page) {
  
  ss_id_regex <- regex(
    r"{
    Shallow\ssample\sID
    \s{15,18}             # Variable whitespace
    (WS\d{6})             # Group 1: worksheet
    _
    (\d{8})               # Group 2: 8 digit DLMS number
    ([^\s]{0,17})         # Group 3: Additional sample info only present for some samples
    # Select any character after the underscore that isn't a space
    # String may be 17 characters. Example: b_0.5_downsampled
    }",
    comments = TRUE
  )
  
  worksheet <- str_extract(page, ss_id_regex, group = 1)
  
  dlms_dna_number <- parse_number(str_extract(page, ss_id_regex, group = 2))
  
  modifier <- str_extract(page, ss_id_regex, group = 3)
  
  sample_id <- str_c(dlms_dna_number, modifier)
  
  shallow_sample_id <- str_c(worksheet, "_", dlms_dna_number, modifier)
  
  return(list(worksheet, dlms_dna_number, modifier, sample_id, shallow_sample_id))
  
}


read_seqone_report <- function(file) {
  seqone_report_text <- pdftools::pdf_text(pdf = file)

  page1 <- seqone_report_text[[1]]

  # HRD score
  
  hrd_score <- get_hrd_score(page = page1, version = "1.1")
  
  assert_that(is.na(hrd_score) == FALSE,
              msg = str_c("Seqone HRD score is NA. File: ", basename(file))
  )
  
  assert_that(hrd_score >= 0, hrd_score <= 1,
              msg = str_c("SeqOne HRD score outside 0-1 range. File: ", basename(file)))

  # HRD status
  
  seqone_hrd_status <- get_hrd_status(page1)

  # LGA
  
  lga <- get_lga(page1)
  
  # LPC
  
  lpc <- get_lpc(page1)

  # CCNE1
  
  ccne1_regex <- regex(
    r"[
    CCNE1\sAmplification
    \s{29}                       # Variable whitespace
    ((\d{1}\.\d{1,2})|(\d{1}))   # Variable number format
    ]",
    comments = TRUE
  )
  
  ccne1 <- parse_number(str_extract(page1, ccne1_regex, group = 1), 
                        locale = locale(decimal_mark = "."))
  
  # RAD51B
  
  rad51b_regex <- regex(
    r"[
    RAD51B\sAmplification
    \s{28}
    ((\d{1}\.\d{1,2})|(\d{1}))
    ]",
    comments = TRUE
  )
  
  rad51b <- parse_number(str_extract(page1, rad51b_regex, group =1),
               locale = locale(decimal_mark = "."))

  # Neoplastic cell content

  seqone_ncc <- get_ncc(page1)
  
  # Coverage
  
  coverage <- get_coverage(page1)

  # Percent mapping
  
  percent_mapping <- get_percent_mapping(page1)

  # Shallow sample ID
  
  id_string <- get_shallow_sample_id(page1)
  
  worksheet <- id_string[[1]]
  
  dlms_dna_number <- id_string[[2]]
  
  modifier <- id_string[[3]]
  
  sample_id <- id_string[[4]]
  
  shallow_sample_id <- id_string[[5]]
  
  # Date
  
  date_regex <- regex(
    r"[
    Date
    \s{32}        # Variable whitespace
    (\w{3,9}      # Month as text. Shortest is May (3), longest September (9)
    \s
    \d{1,2}       # Day (1-31)
    ,\s
    \d{4})        # Year 
    ]",
    comments = TRUE
  )
  
  date <- mdy(str_extract(page1, date_regex, group = 1))
  
  # User
  
  user_regex <- regex(
    r"[
    User
    \s{30,32}         # Variable whitespace
    ([^\s]{10,26})    # Username - any character that isn't space
    ]",
    comments = TRUE
  )
  
  user <- str_extract(page1, user_regex, group = 1)

  # Output table
  
  output <- data.frame(
    "shallow_sample_id" = shallow_sample_id,
    "worksheet" = worksheet,
    "sample_id" = sample_id,
    "dlms_dna_number" = dlms_dna_number,
    "seqone_hrd_score" = hrd_score,
    "seqone_hrd_status" = seqone_hrd_status,
    "lga" = lga,
    "lpc" = lpc,
    "ccne1" = ccne1,
    "rad51b" = rad51b,
    "seqone_ncc" = seqone_ncc,
    "coverage" = coverage,
    "percent_mapping" = percent_mapping,
    "date" = date,
    "user" = user,
    "filename" = basename(file)
  )

  #check_na(output)

  return(output)
}


# Database functions ----------------------------------------------------------------

get_sample_data <- function(sample_vector) {
  
  stopifnot(length(sample_vector) >=1)
  
  sample_query <- paste0("SELECT * FROM MolecularDB.dbo.Samples WHERE LABNO IN (",
                         paste(sample_vector, collapse = ", "),
                         ")")
  
  output_data <- sqlQuery(channel = moldb_connection,
                          query = sample_query) |> 
    janitor::clean_names()
  
  return(output_data)
  
}

# Plot functions --------------------------------------------------------------------

safe_colorblind_palette <- c(
  "#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499",
  "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888"
)

safe_blue <- "#88CCEE"
safe_red <- "#CC6677"
safe_grey <- "#888888"

plot_qc <- function(df = repeat_results, x_var = surname, yvar, outcome) {
  
  sub_df <- df |> 
    filter(shallow_sample_id == "WS133557_21003549")
  
  #title_label <- rlang::englue("{{ yvar }}")
  
  ggplot(df, aes(x = as.character({{ x_var }}), y = {{ yvar }})) +
    geom_point(size = 3, alpha = 0.6,
               aes(colour = {{ outcome }})) +
    scale_colour_manual(name = "",
                        values = c(safe_blue, safe_red, safe_grey)) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          panel.grid = element_blank(),
          legend.position = "bottom") +
    geom_point(
      data = sub_df,
      aes(x = {{ x_var }}, y = {{ yvar }}),
      pch = 21, fill = NA, size = 5, colour = safe_red, stroke = 3
    )
  
}

# Line delineating positive from negative SeqOne results
line_df <- data.frame(
  x =    c(18, 17, 16, 15, 14, 13, 18, 17, 16, 15, 14),
  y =    c(0,   4,  9, 13, 18, 23,  4,  9, 13, 18, 23),
  xend = c(18, 17, 16, 15, 14, 13, 17, 16, 15, 14, 13),
  yend = c(4,   9, 13, 18, 23, 36,  4,  9, 13, 18, 23)
)

plot_lpc_lga <- function(df) {
  
  ggplot(df, aes(x = lga_amended, y = lpc_amended)) +
    geom_point(aes(colour = seqone_hrd_status_amended),
               size = 2, alpha = 0.6) +
    scale_colour_manual(name = "SeqOne HRD Status",
                        values = c(safe_blue, safe_red, safe_grey)) +
    geom_segment(
      data = line_df,
      mapping = aes(x = x, y = y, xend = xend, yend = yend)
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
  
}

plot_variation <- function(df = repeat_variation, yvar) {
  
  ggplot(df, aes(x = , y = {{ yvar }})) +
    geom_boxplot() +
    theme_bw() +
    theme(axis.text.x = element_blank()) +
    ylim(0,30)
  
}

save_hrd_plot <- function(input_plot, input_width = 15, input_height = 12, dpi = 300) {
  # Default inputs allow for presenting a plot as half an A4 page

  ggsave(
    filename = paste0(
      format(Sys.time(), "%Y_%m_%d_%H_%M_%S"),
      "_",
      deparse(substitute(input_plot)), ".png"
    ),
    plot = input_plot,
    device = "png",
    path = hrd_plot_path,
    units = "cm",
    width = input_width,
    height = input_height,
    dpi = 300
  )
}

make_individual_plot <- function(input_sample) {
  output_plot <- compare_results |>
    filter(dlms_dna_number == input_sample) |>
    ggplot(aes(x = worksheet, y = seqone_hrd_score)) +
    geom_point(
      size = 4, alpha = 0.5,
      aes(shape = seqone_hrd_status)
    ) +
    facet_wrap(~dlms_dna_number) +
    theme_bw() +
    labs(
      title = "",
      x = "",
      y = "SeqOne HRD score"
    ) +
    geom_hline(yintercept = 0.50, linetype = "dashed") +
    ylim(0, 1)

  return(output_plot)
}

circle_individual_point <- function(dna_input) {
  results <- results_for_path_block_plot

  output_plot <- path_block_plot +
    geom_point(
      data = results[results$dlms_dna_number == dna_input, ],
      aes(myriad_gi_score, seqone_hrd_score),
      pch = 21, fill = NA, size = 5, colour = safe_red, stroke = 3
    )

  return(output_plot)
}

# Table functions -------------------------------------------------------------------

get_sample_summary_info <- function(input_dna_no) {
  output <- compare_results |>
    filter(dlms_dna_number == input_dna_no) |>
    select(
      worksheet, dlms_dna_number, seqone_hrd_score,
      seqone_hrd_status, lga, lpc, ccne1, rad51b,
      coverage, percent_mapping, million_reads, read_length,
      insert_size, percent_q30, percent_aligned, percent_dups,
      myriad_gi_score, myriad_hrd_status
    )

  return(output)
}

extract_kapa_data <- function(worksheet_number, worksheet_length) {
  description <- grep(
    pattern = worksheet_number,
    x = hs2_library_prep$plate_position,
    value = TRUE
  )

  row_start <- match(description, hs2_library_prep$plate_position) + 1

  row_end <- (row_start + worksheet_length) - 1

  output <- hs2_library_prep[row_start:row_end, ] |>
    select(-starts_with("x")) |>
    mutate(
      worksheet = paste0("WS", worksheet_number),
      shallow_sample_id = paste0(worksheet, "_", lab_number)
    )

  return(output)
}


# Test metric functions -------------------------------------------------------------

count_category <- function(input_table, outcome_column, match_text) {
  
  nrow(input_table |> 
         filter( {{outcome_column }}  == {{ match_text }}))
  
}

compare_tests <- function(input_table, outcome_column) {
  
  true_positives <- count_category(input_table, {{ outcome_column}} , true_pos_text)
  
  true_negatives <- count_category(input_table, {{ outcome_column}} , true_neg_text)
  
  false_positives <- count_category(input_table, {{ outcome_column}} , false_pos_text)
  
  false_negatives <- count_category(input_table, {{ outcome_column}} , false_neg_text)

  incon_positives <- count_category(input_table, {{ outcome_column}} , incon_pos_text)
  
  incon_negatives <- count_category(input_table, {{ outcome_column}} , incon_neg_text)
  
  opa = round(((true_positives + true_negatives) / sum(true_positives, false_positives,
                                                     true_negatives, false_negatives)) * 100, 1)
  
  all_samples = sum(true_positives, true_negatives, false_positives, false_negatives,
                    incon_positives, incon_negatives)
  
  sensitivity = round((true_positives / sum(true_positives, false_positives)) * 100, 1)
  
  specificity = round((true_negatives / sum(true_negatives, false_negatives)) * 100, 1)
  
  inconclusive_rate  = round(((incon_positives + incon_negatives) / all_samples) * 100, 1)
  
  dna_inputs <- nrow(input_table)
  
  unique_samples <- length(unique(input_table$dlms_dna_number))
  
  check <- ifelse(dna_inputs == sum(true_positives, true_negatives, false_positives,
                                    false_negatives, incon_positives, incon_negatives),
                  "DNA inputs match table sum",
                  "Check required")
  
  metrics <- tribble(
    ~"OPA (%)", ~"Sensitivity (%)", ~"Specificity (%)", ~"Unique samples", ~"DNA inputs", ~"Inconclusive rate (%)",
    opa,        sensitivity,       specificity,         unique_samples,    dna_inputs, inconclusive_rate
    )
  
  confusion_matrix <- tribble(
  ~"x",      ~"x",   ~"HRD",          ~"HRP",           ~"Inconclusive",
  "Myriad", "HRD",   true_positives,  false_negatives,  incon_positives,
  "x",      "HRP",   false_positives, true_negatives,   incon_negatives
  )
  
  return(list(confusion_matrix, metrics, check))
  
}

add_version <- function(input_table, version_text) {
  
  input_table |> 
    mutate(Analysis = version_text) |> 
    relocate(Analysis)
  
}
