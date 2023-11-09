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

consistent_text <- "Seqone consistent with Myriad"

inconsistent_text <- "Seqone NOT consistent with Myriad"

inconclusive_text <- "SeqOne inconclusive"

path_block_match_text <- "pathology blocks match" 

path_block_no_match_text <- "pathology blocks DO NOT match"

seqone_status_levels <- c(neg_text, pos_text, incon_text)

consistency_levels <- c(consistent_text, inconsistent_text, inconclusive_text)

# Filepaths -------------------------------------------------------------------------

dev_team_path <- "S:/central shared/Genetics/Mol_Shared/Development.Team/"

seqone_folder <- paste0(dev_team_path, 
                        "SeqOne Homologous Recombination Deficiency Validation/")

hrd_project_path <- paste0(seqone_folder, "HRD R script files/")

hrd_data_path <- paste0(hrd_project_path, "data/")

hrd_output_path <- paste0(hrd_project_path, "outputs/")

hrd_plot_path <- paste0(hrd_project_path, "plots/")

myriad_reports_location <- paste0(hrd_data_path, "myriad_reports/")

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

get_rnumber <- function(page) {
  
  r_number_regex <- regex(
    r"[
    Patient\sID:        # Escape space with backslash s
    \s+                 # 6 spaces
    (R\d{2}-\w{4})      # R number grouped
    ]",
    comments = TRUE
  )
  
  myriad_r_number <- str_extract(page, r_number_regex, group = 1)
  
  assert_that(!is.na(myriad_r_number))
  
  return(myriad_r_number)
}

get_nhs_number <- function(page) {
  
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

  nhs_no_char <- str_extract(page, nhs_no_regex, group = 1)
  
  nhs_no_double <- parse_number(nhs_no_char, locale = locale(grouping_mark = " "))

  assert_that(!is.na(nhs_no_double))
  
  return(nhs_no_double)
  
}

get_name <- function(page) {
  
  patient_name_regex <- regex(
    r"[
    Patient\sName:
    \s+                 # Variable whitespace
    (\D{5,25})          # Names as non-digits with variable lengths
    \n                  # New line marks end of the name
    ]",
    comments = TRUE
  )

  myriad_patient_name <- str_extract(page, patient_name_regex, group = 1)

  assert_that(!is.na(myriad_patient_name))
  
  return(myriad_patient_name)
  
}

get_dob <- function(page) {
  
  dob_regex <- regex(
    r"[
    Date\sof\sBirth:
    \s+                       # Variable whitespace
    (\d{2}/\d{2}/\d{4})       # DOB in dd/mm/yyyy format
    \n                        # Mark end of DOB with newline
    ]",
    comments = TRUE
  )
  
  myriad_dob <- str_extract(page, dob_regex, group = 1)
  
  assert_that(!is.na(myriad_dob))
  
  return(myriad_dob)
  
}

get_path_no <- function(page = page1) {
  
  # Function for extracting pathology number printed on page 1 of Myriad reports
  
  path_block_no_regex <- regex(
    r"[
    Pathology\sNo:
    \s+               # Variable whitespace
    (.{5,25})         # Path blocks can have various formats from different labs
    \n                # Specify end of block number with new line
    ]",
    comments = TRUE
  )
  
  myriad_pathology_number <- str_extract(page, path_block_no_regex, group = 1)
  
  assert_that(!is.na(myriad_pathology_number))
  
  return(myriad_pathology_number)
  
}

get_path_block <- function(page) {
  
  # Function for extracting "pathology block ID" from page 2 of Myriad reports
  
  path_block_regex <- regex(
    r"[
    (Block\(s\)|Specimen\(s\))  # Some say Block(s) whilst others say Specimen(s)
    \sAnalyzed:\s
    (.{5,20})                   # Path blocks can have various formats from different labs
    \n                          # Specify end of block number with new line
    ]",
    comments = TRUE
  )
  
  myriad_pathology_block <- str_extract(page, path_block_regex, group = 2)
  
  assert_that(!is.na(myriad_pathology_block))
  
  return(myriad_pathology_block)
  
}

get_gi_score <- function(page) {
  
  gi_score_regex <- regex(
    r"[
    Patient\sGenomic\sInstability\sScore:
    \s
    (\d{1,3})                               # Score can be 1-3 digits - 1, 10, 100
    \n                                      # Mark end of score with newline
    ]",
    comments = TRUE
  )
  
  gi_score_char <- str_extract(page, gi_score_regex, group = 1)
  
  gi_score_double <- parse_number(gi_score_char)
  
  assert_that(is.na(gi_score_double) == FALSE,
              msg = "GI score is NA")
  
  assert_that(gi_score_double >= 0, gi_score_double <= 100, 
              msg = "GI score outside 0-100 range ")
  
  return(gi_score_double)

}

get_myriad_hrd_status <- function(page) {
  
  hrd_status_regex <- regex(
    r"[
    Myriad\sHRD\sStatus:\s
    (NEGATIVE | POSITIVE)
    ]",
    comments = TRUE
  )

  myriad_hrd_status <- str_extract(page, hrd_status_regex, group = 1)
  
  assert_that(!is.na(myriad_hrd_status))
  
  return(myriad_hrd_status)
  
}

get_myriad_brca_status <- function(page) {

  brca_status_regex <- regex(
    r"[
    Tumor\sMutation\sBRCA1/BRCA2\sStatus:\s
    (NEGATIVE | POSITIVE)
    ]",
    comments = TRUE
    )

  myriad_brca_status <- str_extract(page, brca_status_regex, group = 1)
  
  assert_that(!is.na(myriad_brca_status))
  
  return(myriad_brca_status)
  
}

read_myriad_report <- function(file) {
  
  myriad_report_text <- pdftools::pdf_text(pdf = file)

  page1 <- myriad_report_text[[1]]

  page2 <- myriad_report_text[[2]]

  page3 <- myriad_report_text[[3]]

  output <- data.frame(
    "myriad_r_number" = get_rnumber(page2),
    "myriad_patient_name" = get_name(page1),
    "myriad_dob" = get_dob(page1),
    "nhs_number" = get_nhs_number(page1),
    "myriad_pathology_block_pg1" = get_path_no(page1),
    "myriad_pathology_block_pg2" = get_path_block(page2),
    "myriad_gi_score" = get_gi_score(str_c(page2, page3)),
    "myriad_hrd_status" = get_myriad_hrd_status(page2),
    "myriad_brca_status" = get_myriad_brca_status(page2),
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
    \s+                           # Variable whitespace
    (\d{1}\.\d{1,2} | \d{1})      # Variable formats: 0.99, 0.9, 1
                                  # Use backslash . to specify decimal point
    ]",
    comments = TRUE
  )

  hrd_score_regex_1_2 <- regex(
    r"[
    HRD\sSummary\n
    \s+
    (\d{1}\.\d{1,2} | \d{1})
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

  assert_that(is.na(hrd_score_double) == FALSE,
              msg = "Seqone HRD score is NA")
  
  assert_that(hrd_score_double >= 0, hrd_score_double <= 1,
              msg = "SeqOne HRD score outside 0-1 range")
  
  return(hrd_score_double)
  
}

get_hrd_status <- function(page) {
  
  hrd_status_regex <- regex(
    r"[
    SeqOne\sHRD\sStatus1\s:\s
    (NEGATIVE | POSITIVE | NON-CONCLUSIVE)
    ]",
    comments = TRUE
  )
  
  seqone_hrd_status <- fct(str_extract(page, hrd_status_regex, group = 1),
                           levels = seqone_status_levels)
  
  assert_that(!is.na(seqone_hrd_status))
  
  return(seqone_hrd_status)
  
}

get_lga <- function(page) {
  
  lga_regex <- regex(
    r"[
    LGA\sStatus
    \s+             # Variable whitespace
    (\d{1,2})       # Group LGA status
    ]",
    comments = TRUE
  )
  
  lga <- parse_number(str_extract(page, lga_regex, group = 1))
  
  assert_that(!is.na(lga))
  
  return(lga)
  
}

get_lpc <- function(page) {
  
  lpc_regex <- regex(
    r"[
    LPC\sStatus
    \s+              # Variable whitespace between versions
    (\d{1,2})
    ]",
    comments = TRUE
  )
  
  lpc <- parse_number(str_extract(page, lpc_regex, group = 1))
  
  assert_that(!is.na(lpc))
  
  return(lpc)
  
}

get_ncc <- function(page) {
  
  ncc_regex <- regex(
    r"[
    %\sof\stumoral\scells
    \s+                     # Variable whitespace between versions
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
    \s+                           # Variable whitespace between versions
    (\d{1}\.\d{1,2} | \d{1})      # Variable format 1.57, 1.5, 1
    X
    ]",
    comments = TRUE
  )
  
  coverage <- parse_number(str_extract(page, coverage_regex, group = 1),
                           locale = locale(decimal_mark = "."))
  
  assert_that(!is.na(coverage))
  
  return(coverage)
  
}

get_percent_mapping <- function(page) {
  
  percent_map_regex <- regex(
    r"[
    %\scorrect\smapping
    \s+                           # Variable whitespace between versions
    (\d{2}\.\d{1} | \d{2,3})      # Format 97.2, 97, 100
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
    \s+                   # Variable whitespace
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

get_date <- function(page) {
  
  date_regex <- regex(
    r"[
    Date
    \s+           # Variable whitespace
    (\w{3,9}      # Month as text. Shortest is May (3), longest September (9)
    \s
    \d{1,2}       # Day (1-31)
    ,\s
    \d{4})        # Year 
    ]",
    comments = TRUE
  )
  
  date <- mdy(str_extract(page, date_regex, group = 1))
  
  assert_that(!is.na(date))
  
  return(date)
  
}

get_user <- function(page) {
  
  user_regex <- regex(
    r"[
    User
    \s+               # Variable whitespace
    ([^\s]{10,26})    # Username - any character that isn't space
    ]",
    comments = TRUE
  )
  
  user <- str_extract(page, user_regex, group = 1)
  
  assert_that(!is.na(user))
  
  return(user)
  
}

get_robustness <- function(page, version) {
  
  check_version(version)
  
  robustness_regex <- regex(
    r"[
    Robustness\sof\sgenomic\sinstability
    \s+
    (\d{1}\.\d{1,2} | \d{1})
    ]",
    comments = TRUE
  )
  
  if (version == "1.1") {
    
    robustness <- "not calculated"
    
  }
  
  if (version == "1.2") {
    
    robustness <- parse_number(str_extract(page, robustness_regex, group = 1),
                               locale = locale(decimal_mark = "."))
    
  }
  
  assert_that(!is.na(robustness))
  
  return(robustness)
  
}

get_low_tumour_fraction <- function(page, version) {
  
  check_version(version)
  
  ltf_regex <- regex(
    r"[
    Low\stumor\sfraction
    \s+
    (NORMAL | WARNING)
    ]",
    comments = TRUE
  )
  
  if (version == "1.1") {
    
    low_tumour_fraction <- "not calculated"
    
  }
  
  if (version == "1.2") {
    
    low_tumour_fraction <- str_extract(page, ltf_regex, group = 1)
    
  }
  
  assert_that(!is.na(low_tumour_fraction))
  
  return(low_tumour_fraction)
  
}

get_ccne1_rad51b <- function(page, version) {
  
  check_version(version)
  
  ccne1_regex_1_1 <- regex(
    r"[
    CCNE1\sAmplification
    \s+                              # Variable whitespace
    (\d{1,2}\.\d{1,2} | \d{1,2})     # Variable number format
    ]",
    comments = TRUE
  )
  
  rad51b_regex_1_1 <- regex(
    r"[
    RAD51B\sAmplification
    \s+
    (\d{1,2}\.\d{1,2} | \d{1,2})
    ]",
    comments = TRUE
  )
  
  ccne1_regex_1_2 <- regex(
    r"[
    RAD51B
    \n\n
    \s+
    (\d{1,2}\.\d{1,2} | \d{1,2})    # CCNE1 regex
    \s+
    (\d{1,2}\.\d{1,2} | \d{1,2})    # RAD51B regex
    ]",
    comments = TRUE
  )
  
  if (version == "1.1") {
    
    ccne1 <- parse_number(str_extract(page, ccne1_regex_1_1, group = 1),
               locale = locale(decimal_mark = "."))
    
    rad51b <-  parse_number(str_extract(page, rad51b_regex_1_1, group = 1),
                                      locale = locale(decimal_mark = "."))
    
  }
  
  if (version == "1.2") {
    
    ccne1 <- parse_number(str_extract(page, ccne1_regex_1_2, group = 1),
                          locale = locale(decimal_mark = "."))
    
    rad51b <- parse_number(str_extract(page, ccne1_regex_1_2, group = 2),
                           locale = locale(decimal_mark = "."))
    
  }
  
  assert_that(!is.na(ccne1))
  
  assert_that(!is.na(rad51b))
  
  return(list(ccne1, rad51b))
  
}

read_seqone_report <- function(file, version) {
  
  check_version(version)
  
  seqone_report_text <- pdftools::pdf_text(pdf = file)

  page1 <- seqone_report_text[[1]]

  # Output table
  
  output <- data.frame(
    "shallow_sample_id" = get_shallow_sample_id(page1)[[5]],
    "worksheet" = get_shallow_sample_id(page1)[[1]],
    "sample_id" = get_shallow_sample_id(page1)[[4]],
    "dlms_dna_number" = get_shallow_sample_id(page1)[[2]],
    "seqone_hrd_score" = get_hrd_score(page = page1, version = version),
    "seqone_hrd_status" = get_hrd_status(page1),
    "lga" = get_lga(page1),
    "lpc" = get_lpc(page1),
    "ccne1" = get_ccne1_rad51b(page1, version)[[1]],
    "rad51b" = get_ccne1_rad51b(page1, version)[[2]],
    "seqone_ncc" = get_ncc(page1),
    "coverage" = get_coverage(page1),
    "percent_mapping" = get_percent_mapping(page1),
    "robustness" = get_robustness(page1, version),
    "low_tumour_fraction" = get_low_tumour_fraction(page1, version),
    "date" = get_date(page1),
    "user" = get_user(page1),
    "filename" = basename(file),
    "version" = version
  )

  check_na(output)

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

safe_blue <- "#88CCEE"
safe_red <- "#CC6677"
safe_grey <- "#888888"

plot_qc <- function(df = filtered_results, x_var = shallow_sample_id, yvar, 
                    outcome = outcome_binary, alpha_number = 0.6) {
  
  ggplot(df, aes(x = as.character({{ x_var }}), y = {{ yvar }})) +
    geom_point(size = 3, alpha = alpha_number,
               aes(colour = {{ outcome }})) +
    scale_colour_manual(name = "",
                        values = c(safe_blue, safe_red, safe_grey,
                                   safe_grey)) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          panel.grid = element_blank(),
          legend.position = "bottom") +
    labs(x = "")
  
}

# Line delineating positive from negative SeqOne results
line_df <- data.frame(
  x =    c(18, 17, 16, 15, 14, 13, 18, 17, 16, 15, 14),
  y =    c(0,   4,  9, 13, 18, 23,  4,  9, 13, 18, 23),
  xend = c(18, 17, 16, 15, 14, 13, 17, 16, 15, 14, 13),
  yend = c(4,   9, 13, 18, 23, 36,  4,  9, 13, 18, 23)
)

plot_lpc_lga <- function(df) {
  
  ggplot(df, aes(x = lga, y = lpc)) +
    geom_point(aes(colour = seqone_hrd_status),
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
    theme(axis.text.x = element_blank()) 
  
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

# Table functions -------------------------------------------------------------------

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

make_robustness_table <- function(filtered_df) {
  
  output <- compare_results |> 
    filter(shallow_sample_id %in% filtered_df$shallow_sample_id) |> 
    select(dlms_dna_number, seqone_hrd_status,
           myriad_hrd_status, robustness,
           path_block_manual_check, version) |> 
    pivot_wider(id_cols = c(dlms_dna_number, myriad_hrd_status,
                            path_block_manual_check),
                names_from = version,
                values_from = c(seqone_hrd_status, robustness)) |> 
    select(dlms_dna_number, robustness_1.2,
           seqone_hrd_status_1.2, seqone_hrd_status_1.1, 
           myriad_hrd_status, 
           path_block_manual_check) |> 
    arrange(path_block_manual_check)
  
  return(output)
  
}

# Test metric functions -------------------------------------------------------------

compare_tests <- function(input_table) {
  
  assert_that("outcome" %in% colnames(input_table))
  
  assert_that(anyNA(input_table$outcome) == FALSE)
  
  assert_that(!"other" %in% input_table$outcome)
  
  true_positives <- nrow(input_table[input_table$outcome == true_pos_text, ])
  
  true_negatives <- nrow(input_table[input_table$outcome == true_neg_text, ])

  false_positives <- nrow(input_table[input_table$outcome == false_pos_text, ])
  
  false_negatives <- nrow(input_table[input_table$outcome == false_neg_text, ])

  incon_positives <- nrow(input_table[input_table$outcome == incon_pos_text, ])
  
  incon_negatives <- nrow(input_table[input_table$outcome == incon_neg_text, ])
  
  opa = round(((true_positives + true_negatives) / sum(true_positives, false_positives,
                                                     true_negatives, false_negatives)) * 100, 1)
  
  all_samples = sum(true_positives, true_negatives, false_positives, false_negatives,
                    incon_positives, incon_negatives)
  
  sensitivity = round((true_positives / sum(true_positives, false_negatives)) * 100, 1)
  
  specificity = round((true_negatives / sum(true_negatives, false_positives)) * 100, 1)
  
  inconclusive_rate  = round(((incon_positives + incon_negatives) / all_samples) * 100, 1)
  
  dna_inputs <- nrow(input_table)
  
  unique_samples <- length(unique(input_table$dlms_dna_number))
  
  check <- ifelse(dna_inputs == sum(true_positives, true_negatives, false_positives,
                                    false_negatives, incon_positives, incon_negatives),
                  "DNA inputs match table sum",
                  "Check required")
  
  metrics <- tribble(
    ~"OPA (%)", ~"Sensitivity (%)", ~"Specificity (%)", ~"Unique samples", ~"DNA inputs", ~"Inconclusive rate (%)",
    opa,        sensitivity,       specificity,         unique_samples,    dna_inputs,     inconclusive_rate
    )
  
  confusion_matrix <- tribble(
  ~"x",      ~"x",       ~"HRD (+)",      ~"HRP (-)",        ~"Inconclusive",
  "Myriad", "HRD (+)",   true_positives,  false_negatives,  incon_positives,
  "x",      "HRP (-)",   false_positives, true_negatives,   incon_negatives
  )
  
  return(list(confusion_matrix, metrics, check))
  
}

add_version <- function(input_table, version_text) {
  
  input_table |> 
    mutate(Analysis = version_text) |> 
    relocate(Analysis)
  
}
