# Collate SeqOne Audit Data

# Libraries and packages ------------------------------------------------------------

library("tidyverse")
library("here")
library("readxl")

source(here("functions/hrd_functions.R"))

data_config <- config::get("data_filepath")

# Collate validation PDF data -------------------------------------------------------

validation_reports <- list.files(path = paste0(data_config, 
                                 "validation/DOC6192_seqone/seqone_reports_v1_2/"),
                                 full.names = TRUE,
                                 pattern = "*.pdf")

collated_validation_pdf_data <- validation_reports |>
  map(\(validation_reports) read_seqone_report(
    file = validation_reports,
    version = "1.2"
  )) |>
  list_rbind()

if(anyNA.data.frame(collated_validation_pdf_data) == TRUE |
   nrow(collated_validation_pdf_data) == 0) {
  stop()
}

write.csv(collated_validation_pdf_data, 
          file = paste0(data_config, 
                 "validation/DOC6192_seqone/collated_data/", 
                 "collated_validation_pdf_data.csv"),
          row.names = FALSE)

# Collate QIAsymphony validation samples ---------------------------------------------------------

# Worksheet WS136827 was worksheet of Qiasymphony validation samples 
# that was run after the initial validation ended and before the live service started.

qiasymphony_reports <- list.files(paste0(data_config, 
                                      "validation/DOC6255_qiasymphony/seqone_reports/"),
                               full.names = TRUE,
                               pattern = "*.pdf")

qiasymphony_pdf_data <- qiasymphony_reports |>
  map(\(qiasymphony_reports) read_seqone_report(
    file = qiasymphony_reports,
    version = "1.2"
  )) |>
  list_rbind()

if(anyNA.data.frame(qiasymphony_pdf_data) == TRUE |
   nrow(qiasymphony_pdf_data) == 0) {
  stop()
}

write.csv(qiasymphony_pdf_data, 
          file = paste0(data_config, 
                        "validation/DOC6255_qiasymphony/collated_data/", 
                        "qiasymphony_pdf_data.csv"),
          row.names = FALSE)

# Load INC9096 data -----------------------------------------------------------------

# This is the data for the samples that had abnormally low confidence
# in genomic instability scores due to a pipeline error (see incident INC9096).

inc9096_pdfs <- list.files(paste0(data_config, 
                                   "live_service/INC9096/data/"),
                            full.names = TRUE,
                            pattern = ".pdf")

inc9096_pdf_data <- inc9096_pdfs |> 
  map(\(inc9096_pdfs) read_seqone_report(file = inc9096_pdfs,
                                          version = "1.2")) |> 
  list_rbind() |> 
  mutate(worksheet = "WS140359_1") |> 
  filter(date == "2024-04-01")

if(anyNA.data.frame(inc9096_pdf_data) == TRUE |
   nrow(inc9096_pdf_data) == 0) {
  stop()
}

write.csv(inc9096_pdf_data, 
          paste0(data_config, 
                 "live_service/INC9096/collated_data/inc9096_pdf_data.csv"),
          row.names = FALSE)

# Collate live service PDF data -----------------------------------------------------

live_ws_excel <- read_excel(path = paste0(data_config, 
                                          "live_service/live_service_worksheets.xlsx"))

live_ws <- list(live_ws_excel$worksheet)

live_ws_filepaths_pdf <- live_ws |> 
  map(\(live_ws) find_hrd_files(live_ws, filetype = ".pdf")) |> 
  flatten()

file.copy(from = live_ws_filepaths_pdf,
          to = paste0(data_config, 
                      "live_service/service/pdf_reports/"))

data_folder_live_pdfs <- list.files(path = paste0(data_config, 
                                                  "live_service/service/pdf_reports/"),
                                         pattern = "*.pdf",
                                         full.names = TRUE)

collated_live_pdf_data <- data_folder_live_pdfs |>
  map(\(data_folder_live_pdfs) read_seqone_report(
    file = data_folder_live_pdfs,
    version = "1.2"
  )) |>
  list_rbind() |> 
  mutate(worksheet = ifelse(worksheet == "WS140359", "WS140359_2", worksheet))

if(nrow(collated_live_pdf_data) == 0) {
  stop()
}

write.csv(collated_live_pdf_data, 
          paste0(data_config, 
                 "live_service/service/collated_data/",
                 "collated_live_pdf_data.csv"),
          row.names = FALSE)

# Collate live service CSV data -----------------------------------------------------

live_ws_filepaths_csv <- live_ws |> 
  map(\(live_ws) find_hrd_files(live_ws, filetype = ".csv")) |> 
  flatten()

file.copy(from = live_ws_filepaths_csv,
          to = paste0(data_config, "live_service/service/",
                      "csv_reports"))

data_folder_live_csv <- list.files(path = paste0(data_config, 
                                                 "live_service/service/",
                                               "csv_reports"),
                                          pattern = "hrd-results.*csv",
                                          full.names = TRUE)

collated_live_csv_data <- data_folder_live_csv |> 
  map(\(data_folder_live_csv) read_seqone_csv(file = data_folder_live_csv)) |> 
  list_rbind() |> 
  mutate(worksheet = str_extract(string = sample,
                                 pattern = "(WS[0-9]{6})_(\\d{8})",
                                 group = 1),
         labno = str_extract(string = sample,
                             pattern = "(WS[0-9]{6})_(\\d{8})",
                             group = 2))

if(nrow(collated_live_csv_data) == 0) {
  stop()
}

write.csv(collated_live_csv_data, 
          paste0(data_config, 
                 "live_service/service/collated_data/",
                 "collated_live_csv_data.csv"),
          row.names = FALSE)
