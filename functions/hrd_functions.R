################################################################################
# HRD Project Functions
################################################################################

source("scripts/hrd_filepaths.R")

library(pdftools)
library(tidyverse)
library(assertthat)

##################################################
# CSV Timestamp
##################################################

export_timestamp <- function(filepath, input) {
  
  write.csv(input, 
            file = paste0(filepath,
                          format(Sys.time(), "%Y_%m_%d_%H_%M_%S"),
                          "_",
                          deparse(substitute(input)), ".csv"),
            row.names = FALSE)
}

##################################################
# Plot Functions
##################################################

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

safe_blue <- "#88CCEE"
safe_red <- "#CC6677"

save_hrd_plot <- function(input_plot, input_width = 15, input_height = 12, dpi = 300) {
  
  # Default inputs allow for presenting a plot as half an A4 page
  
  ggsave(filename = paste0(format(Sys.time(), "%Y_%m_%d_%H_%M_%S"),
                           "_",
                           deparse(substitute(input_plot)), ".png"),
         plot = input_plot,
         device = "png",
         path = hrd_plot_path,
         units = "cm",
         width = input_width,
         height = input_height,
         dpi = 300)
  
}

make_individual_plot <- function(input_sample) {
  
  output_plot <- compare_results |> 
    filter(dlms_dna_number == input_sample) |> 
    ggplot(aes(x = worksheet, y = seqone_hrd_score)) +
    geom_point(size = 4, alpha = 0.5, 
               aes(shape = seqone_hrd_status)) +
    facet_wrap(~dlms_dna_number) +
    theme_bw() +
    labs(title = "",
         x = "",
         y = "SeqOne HRD score") +
    geom_hline(yintercept = 0.50, linetype = "dashed") +
    ylim(0, 1)
  
  return(output_plot)
  
}

circle_individual_point <- function(dna_input) {
  
  output_plot <- path_block_plot +
    geom_point(data=results_for_path_block_plot[results_for_path_block_plot$dlms_dna_number == dna_input,], 
               aes(myriad_gi_score, seqone_hrd_score),
               pch=21, fill=NA, size=5, colour = safe_red, stroke=3)
  
  return(output_plot)
  
}

##################################################
# Table Functions
##################################################

get_sample_summary_info <- function(input_dna_no) {
  
  output <- compare_results |> 
    filter(dlms_dna_number == input_dna_no)  |> 
    select(worksheet, dlms_dna_number, seqone_hrd_score,
           seqone_hrd_status, lga, lpc, ccne1, rad51b, 
           coverage.x, percent_mapping, million_reads, read_length, 
           insert_size, percent_q30, percent_aligned, percent_dups,
           myriad_gi_score, myriad_hrd_status)
  
  return(output)
  
}

extract_kapa_data <- function(worksheet_number, worksheet_length) {
  
  description <- grep(pattern = worksheet_number,
                      x = hs2_library_prep$plate_position,
                      value = TRUE)
  
  row_start <- match(description, hs2_library_prep$plate_position) + 1 
  
  row_end <- (row_start + worksheet_length) - 1
  
  output <- hs2_library_prep[row_start:row_end,] |>
    select(-starts_with("x")) |>
    mutate(worksheet = paste0("WS", worksheet_number),
           shallow_sample_id = paste0(worksheet, "_", lab_number))
  
  return(output)
  
}

##################################################
# PDF Functions
##################################################

check_na <- function(input_table) {
  
  assert_that(sum(is.na(input_table)) == 0, msg = "NA values present")
  
}


read_myriad_report <- function(filepath, file) {
  
  # Use pdf_text to read PDF as a single string per page
  
  myriad_report_text <- pdftools::pdf_text(pdf = paste0(filepath, file))
  
  page1 <- myriad_report_text[[1]]
  
  page2 <- myriad_report_text[[2]]
  
  page3 <- myriad_report_text[[3]]
  
  # iGene R Number
  
  myriad_r_number <- sub(x = page2,
                         pattern = ".+Patient ID:.{6}(R\\d{2}-\\w{4}).+",
                         replacement = "\\1")
  
  # NHS Number
  # Note: size of whitespace between "No:" and number can vary
  
  nhs_number <- sub(x = page1,
                    pattern = ".+NHS No:.{10,20}(\\d{3}.{1}\\d{3}.{1}\\d{4}).+",
                    replacement = "\\1")
  
  nhs_number_mod <- as.numeric(gsub(pattern = "(\\D)", "",
                                    nhs_number))
  
  
  # Patient name
  myriad_patient_name <- sub(x = page1,
                             pattern = ".+Patient Name:\\s{10,15}(\\D{5,25})\n.+",
                             replacement = "\\1")
  
  # Patient date of birth
  myriad_dob <- sub(x = page1,
                    pattern = ".+Date of Birth:\\s{10,15}(\\d{2}/\\d{2}/\\d{4})\n.+",
                    replacement = "\\1")
  
  # Pathology Block
  # Note: some reports say "Block(s)" whilst others say "Specimen(s)"
  
  myriad_pathology_block_pg2 <- gsub("[\n]", "", sub(x = page2,
                                                     pattern = ".+(Block\\(s\\)|Specimen\\(s\\)) Analyzed:(.{5,20})\n\n.+",
                                                     replacement = "\\2"))
  
  myriad_pathology_block_pg1 <- gsub("[\n]", "", sub(x = page1,
                                                     pattern = ".+Pathology No:\\s{10,30}(.{5,25})(\\n\\n\\n\\n\\n|\\n).+",
                                                     replacement = "\\1"))
  
  # Genomic Instability Score
  # Note: score can be 1 digit (i.e. "7") or 2 ("73")
  
  gi_score_regex <- ".+Patient Genomic Instability Score: (\\d{1,2}).+"
  
  gi_score_pg2 <- as.numeric(sub(x = page2,
                                 pattern = gi_score_regex,
                                 replacement = "\\1"))
  
  # Some reports have GI score on the third page
  
  gi_score_pg3 <- as.numeric(sub(x = page3,
                                 pattern = gi_score_regex,
                                 replacement = "\\1"))
  
  # Pick correct GI score
  
  myriad_gi_score <- ifelse(str_length(gi_score_pg2) %in% c(1,2),
                            gi_score_pg2,
                            ifelse(str_length(gi_score_pg3) %in% c(1,2),
                                   gi_score_pg3, "NULL"))
  
  assert_that(myriad_gi_score >= 0, myriad_gi_score <= 100, msg = "GI score outside 0-100 range")
  
  
  myriad_hrd_status <- sub(x = page2,
                           pattern = ".+Myriad HRD Status:.(\\D{8}).+",
                           replacement = "\\1")
  
  assert_that(myriad_hrd_status %in% c("POSITIVE", "NEGATIVE"))
  
  myriad_brca_status <- sub(x = page2,
                            pattern = ".+Tumor Mutation BRCA1/BRCA2 Status:.(\\D{8}).+",
                            replacement = "\\1")
  
  
  output <- data.frame("myriad_r_number" = myriad_r_number,
                       "myriad_patient_name" = myriad_patient_name,
                       "myriad_dob" = myriad_dob,
                       "nhs_number" = nhs_number_mod,
                       "myriad_pathology_block_pg1" = myriad_pathology_block_pg1,
                       "myriad_pathology_block_pg2" = myriad_pathology_block_pg2,
                       "myriad_gi_score" = myriad_gi_score,
                       "myriad_hrd_status" = myriad_hrd_status,
                       "myriad_brca_status" = myriad_brca_status,
                       "myriad_filename" = file)
  
  check_na(output)
  
  return(output)
  
}

grep_seqone_text <- function(input_regex, page) {
  
  output <- sub(pattern = input_regex,
                x = page,
                replacement = "\\1")
  
  return(output)
  
}

read_seqone_report <- function(filepath, file) {
  
  seqone_report_text <- pdftools::pdf_text(pdf = paste0(filepath, file))
  
  page_1 <- seqone_report_text[[1]]
  
  seqone_hrd_score <- as.numeric(grep_seqone_text(".+CLASS\n\\s{78,83}((0.\\d{2})|(0.\\d{1})|(\\d{1})).+", 
                                                  page_1))
  
  assert_that(seqone_hrd_score >= 0, seqone_hrd_score <=1, msg = "SeqOne HRD score outside 0-1 range")
  
  seqone_hrd_status <- grep_seqone_text(".+SeqOne HRD Status1 : (\\D{8}).+", page_1)
  
  assert_that(seqone_hrd_status %in% c("POSITIVE", "NEGATIVE"), msg = "SeqOne HRD score not dichotomous")
  
  lga <- as.numeric(grep_seqone_text(".+LGA Status.{38}(\\d{1,2}).+", page_1))
  
  lpc <- as.numeric(grep_seqone_text(".+LPC Status.{38}(\\d{1,2}).+", page_1))
  
  ccne1 <- as.numeric(grep_seqone_text(".+CCNE1 Amplification.{29}((\\d{1}.\\d{2})|\\d{1}).+", page_1))
  
  rad51b <- as.numeric(grep_seqone_text(".+RAD51B Amplification.{28}(\\d{1}.\\d{1,2}|\\d{1}).+", page_1))
  
  seqone_ncc <- as.numeric(grep_seqone_text(".+% of tumoral cells.{23,24}(\\d{2})%.+", page_1))
  
  assert_that(seqone_ncc >= 20, seqone_ncc <= 100)
  
  coverage <- as.numeric(grep_seqone_text(".+Coverage\\s{33,34}(.{1,4})X.+", page_1))
  
  percent_mapping <- as.numeric(grep_seqone_text(".+% correct mapping.{24,25}((\\d{2}.\\d{1})|(\\d{2}))%.+", page_1))
  
  assert_that(percent_mapping >= 0, percent_mapping <= 100)
  
  shallow_sample_id <- trimws(grep_seqone_text(".+Shallow sample ID\\s{15,18}(WS\\d{6}_.{8,26}).+", page_1),
                              which = "right")
  
  sample_id <- sub(pattern = "WS\\d{6}_(.{8,26})",
                   x = shallow_sample_id,
                   replacement = "\\1")
  
  worksheet <- sub(pattern = "^(WS\\d{6}).+",
                   x = shallow_sample_id,
                   replacement = "\\1")
  
  dlms_dna_number <- as.numeric(sub(pattern = "^(\\d{8}).+",
                                    replacement = "\\1",
                                    x = sample_id))
  
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
    "rad51b" =rad51b,
    "seqone_ncc" = seqone_ncc,
    "coverage" = coverage,
    "percent_mapping" = percent_mapping,
    "date" = date,
    "user" = user,
    "filename" = file)
  
  check_na(output)
  
  return(output)
  
}

##################################################
# Collation Functions
##################################################

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
  collated_seqone_check <- output  |>
    # Combination of sample ID and run date should be unique
    mutate(sample_date = paste0(shallow_sample_id, " ", date))
  
  stopifnot(!duplicated(collated_seqone_check$sample_date))
  
  rm(collated_seqone_check)
  
  return(output)

}

##################################################