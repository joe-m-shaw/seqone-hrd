# HRD Project Functions

source("scripts/hrd_filepaths.R")

library(pdftools)
library(tidyverse)
library(assertthat)


# CSV timestamp ---------------------------------------------------------------------

export_timestamp <- function(filepath, input) {
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



safe_colorblind_palette <- c(
  "#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499",
  "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888"
)

safe_blue <- "#88CCEE"
safe_red <- "#CC6677"

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
      coverage.x, percent_mapping, million_reads, read_length,
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


# PDF functions ---------------------------------------------------------------------

check_na <- function(input_table) {
  assert_that(sum(is.na(input_table)) == 0,
    msg = paste0("NA values present: file ", file)
  )
}


read_myriad_report <- function(filepath, file) {
  # Use pdf_text to read PDF as a single string per page

  myriad_report_text <- pdftools::pdf_text(pdf = paste0(filepath, file))

  page1 <- myriad_report_text[[1]]

  page2 <- myriad_report_text[[2]]

  page3 <- myriad_report_text[[3]]

  # iGene R Number
  
  r_number_regex <- regex(
    r"[
    Patient\sID:        # Escape space with \s
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
    .{10,20}                    # Size of whitespace between "No:" and number can vary
                                # Use . instead of \s as samples can have a leading
                                # 0 before NHS number (example: R22-031L)
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
    (Block\(s\)|Specimen\(s\))  # Some say "Block(s)" whilst others say "Specimen(s)"
    \sAnalyzed:\s
    (.{5,20})                   # Path blocks can have various formats from different labs
    \n                          # Specify end of block number with new line
    ]",
    comments = TRUE
  )
  
  myriad_pathology_block_pg2 <- str_extract(page2, path_block_2_regex, group = 2)

  gi_score_regex <- regex(
    r"[
    Patient\sGenomic\sInstability\sScore:
    \s
    (\d{1,2})                               # Score can be 1 digit (i.e. "7") or 2 ("73")
    \n                                      # Mark end of score with newline
    ]",
    comments = TRUE
  )
  
  gi_score_char <- str_extract(str_c(page2, page3), gi_score_regex, group = 1)
  
  gi_score_double <- parse_number(gi_score_char)

  assert_that(is.na(gi_score_double) == FALSE,
    msg = paste0("GI score is NA: file ", file)
  )

  assert_that(gi_score_double >= 0, gi_score_double <= 100, msg = paste0(
    "GI score outside 0-100 range: file ", file
  ))

  hrd_status_regex <- regex(
    r"[
    Myriad\sHRD\sStatus:\s
    (NEGATIVE | POSITIVE)
    ]",
    comments = TRUE
  )
  
  myriad_hrd_status <- str_extract(page2, hrd_status_regex, group = 1)

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
    "myriad_filename" = file
  )

  check_na(output)

  return(output)
}

grep_seqone_text <- function(input_regex, page) {
  output <- sub(
    pattern = input_regex,
    x = page,
    replacement = "\\1"
  )

  return(output)
}

read_seqone_report <- function(filepath, file) {
  seqone_report_text <- pdftools::pdf_text(pdf = paste0(filepath, file))

  page_1 <- seqone_report_text[[1]]

  seqone_hrd_score <- as.numeric(grep_seqone_text(
    ".+CLASS\n\\s{78,83}((0.\\d{2})|(0.\\d{1})|(\\d{1})).+",
    page_1
  ))

  assert_that(is.na(seqone_hrd_score) == FALSE,
    msg = paste0("Seqone HRD score is NA: file", file)
  )

  assert_that(seqone_hrd_score >= 0, seqone_hrd_score <= 1,
    msg = "SeqOne HRD score outside 0-1 range"
  )

  seqone_hrd_status <- grep_seqone_text(".+SeqOne HRD Status1 : (\\D{8}).+", page_1)

  assert_that(seqone_hrd_status %in% c("POSITIVE", "NEGATIVE"),
    msg = "SeqOne HRD score not dichotomous"
  )

  lga <- as.numeric(grep_seqone_text(".+LGA Status.{38}(\\d{1,2}).+", page_1))

  lpc <- as.numeric(grep_seqone_text(".+LPC Status.{38}(\\d{1,2}).+", page_1))

  ccne1 <- as.numeric(grep_seqone_text(
    ".+CCNE1 Amplification.{29}((\\d{1}.\\d{2})|\\d{1}).+", page_1
  ))

  rad51b <- as.numeric(grep_seqone_text(
    ".+RAD51B Amplification.{28}(\\d{1}.\\d{1,2}|\\d{1}).+", page_1
  ))

  seqone_ncc <- as.numeric(grep_seqone_text(
    ".+% of tumoral cells.{23,24}(\\d{2})%.+", page_1
  ))

  assert_that(seqone_ncc >= 20, seqone_ncc <= 100)

  coverage <- as.numeric(grep_seqone_text(
    ".+Coverage\\s{33,34}(.{1,4})X.+", page_1
  ))

  percent_mapping <- as.numeric(grep_seqone_text(
    ".+% correct mapping.{24,25}((\\d{2}.\\d{1})|(\\d{2}))%.+", page_1
  ))

  assert_that(percent_mapping >= 0, percent_mapping <= 100)

  shallow_sample_id <- trimws(
    grep_seqone_text(
      ".+Shallow sample ID\\s{15,18}(WS\\d{6}_.{8,26}).+", page_1
    ),
    which = "right"
  )

  sample_id <- sub(
    pattern = "WS\\d{6}_(.{8,26})",
    x = shallow_sample_id,
    replacement = "\\1"
  )

  worksheet <- sub(
    pattern = "^(WS\\d{6}).+",
    x = shallow_sample_id,
    replacement = "\\1"
  )

  dlms_dna_number <- as.numeric(sub(
    pattern = "^(\\d{8}).+",
    replacement = "\\1",
    x = sample_id
  ))

  date <- grep_seqone_text(".+Date.{32}(\\D{3,9}.\\d{1,2}..\\d{4}).+", page_1)

  user <- grep_seqone_text(".+User\\s{30,32}(.{10,26})\\s+.+", page_1)

  output <- data.frame(
    "shallow_sample_id" = shallow_sample_id,
    "worksheet" = worksheet,
    "sample_id" = sample_id,
    "dlms_dna_number" = dlms_dna_number,
    "seqone_hrd_score" = seqone_hrd_score,
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
    "filename" = file
  )

  check_na(output)

  return(output)
}

# Collation functions ---------------------------------------------------------------

collate_myriad_reports <- function() {
  myriad_report_files <- list.files(myriad_reports_location)

  output <- data.frame()

  for (i in myriad_report_files) {
    tmp_output <- read_myriad_report(myriad_reports_location, i)

    output <- rbind(output, tmp_output)

    rm(tmp_output)
  }

  # Check all files included
  stopifnot(setdiff(myriad_report_files, output$myriad_filename) == 0)

  return(output)
}

collate_seqone_reports <- function() {
  seqone_reports <- list.files(seqone_report_location)

  output <- data.frame()

  for (i in seqone_reports) {
    tmp_output <- read_seqone_report(seqone_report_location, i)

    output <- rbind(output, tmp_output)

    rm(tmp_output)
  }

  # Check all files collated
  stopifnot(setdiff(seqone_reports, output$filename) == 0)

  # Check to make sure a report wasn't saved twice
  collated_seqone_check <- output |>
    # Combination of sample ID and run date should be unique
    mutate(sample_date = paste0(shallow_sample_id, " ", date))

  stopifnot(!duplicated(collated_seqone_check$sample_date))

  rm(collated_seqone_check)

  return(output)
}

# Test metric functions -------------------------------------------------------------

perform_sensitivity_calcs <- function(input_table) {
  
  required_columns <- c("seqone_hrd_status_amended", "hrd_status_check_amended",
                        "myriad_hrd_status")
  
  stopifnot(required_columns %in% colnames(input_table))
  
  input_table_mod <- input_table |> 
    mutate(
      
      category = case_when(
        
        seqone_hrd_status_amended == "POSITIVE" &
          hrd_status_check_amended == consistent_text ~"true positive",
        
        seqone_hrd_status_amended == "NEGATIVE" &
          hrd_status_check_amended == consistent_text ~"true negative",
        
        seqone_hrd_status_amended == "POSITIVE" &
          hrd_status_check_amended == inconsistent_text ~"false positive",
        
        seqone_hrd_status_amended == "NEGATIVE" &
          hrd_status_check_amended == inconsistent_text ~"false negative",
        
        seqone_hrd_status_amended == "NON-CONCLUSIVE" &
          myriad_hrd_status == "POSITIVE" ~"inconclusive positive",
        
        seqone_hrd_status_amended == "NON-CONCLUSIVE" &
          myriad_hrd_status == "NEGATIVE" ~"inconclusive negative"
      )
      
    )
  
  counts <- input_table_mod |> 
    count(category)
  
  true_positives <- sum(input_table_mod$category == "true positive")
  
  true_negatives <- sum(input_table_mod$category == "true negative")
  
  false_negatives <- sum(input_table_mod$category == "false negative")
  
  false_positives <- sum(input_table_mod$category == "false positive")
  
  inconclusive_negatives <- sum(input_table_mod$category == "inconclusive negative")
  
  inconclusive_positives <- sum(input_table_mod$category == "inconclusive positive")
  
  overall_stats <- tribble(
    
    ~col, ~seqone_positive, ~seqone_negative, ~seqone_inconclusive,
    "myriad_positive", true_positives, false_negatives, inconclusive_positives,
    "myriad_negative", false_positives, true_negatives, inconclusive_negatives
  )
  
  results_minus_inconclusives <- sum(true_positives, true_negatives, 
                                     false_positives, false_negatives)
  
  opa <- ((true_positives + true_negatives) / results_minus_inconclusives) * 100
  
  sensitivity <- (true_positives / sum(true_positives, false_positives)) * 100
  
  specificity <- (true_negatives / sum(true_negatives, false_negatives)) * 100
  
  samples <- length(unique(input_table_mod$dlms_dna_number))
  
  dna_inputs <- sum(true_positives, true_negatives, 
                    false_positives, false_negatives,
                    inconclusive_negatives, inconclusive_positives)
  
  print(overall_stats)
  
  print(paste0("OPA = ", round(opa, 2), "%"))
  
  print(paste0("Sensitivity = ", round(sensitivity, 2), "%"))
  
  print(paste0("Specificity = ", round(specificity, 2), "%"))
  
  print(paste0("DNA inputs: ", dna_inputs))
  
  print(paste0("Unique samples = ", samples))
  
}

