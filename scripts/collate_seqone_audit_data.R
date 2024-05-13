# Collate SeqOne Audit Data

# Libraries and packages ------------------------------------------------------------

library("tidyverse")
library("here")
library("readxl")

source(here::here("functions/hrd_functions.R"))

# Collate validation PDF data -------------------------------------------------------

validation_reports <- list.files(here::here("data/seqone_reports_v1_2/"),
                                 full.names = TRUE,
                                 pattern = "*.pdf")

collated_validation_pdf_data <- validation_reports |>
  map(\(validation_reports) read_seqone_report(
    file = validation_reports,
    version = "1.2"
  )) |>
  list_rbind()

write.csv(collated_validation_pdf_data, 
          here::here("data/seqone_collated_audit_data/collated_validation_pdf_data.csv"),
          row.names = FALSE)

# Collate WS136827 ------------------------------------------------------------------

# Worksheet WS136827 was worksheet of Qiasymphony validation samples 
# that was run after the validation ended and before the live service started.

ws136827_reports <- list.files(here::here("data/seq_one_reports_WS136827/"),
                               full.names = TRUE,
                               pattern = "*.pdf")

ws136827_pdf_data <- ws136827_reports |>
  map(\(ws136827_reports) read_seqone_report(
    file = ws136827_reports,
    version = "1.2"
  )) |>
  list_rbind()

write.csv(ws136827_pdf_data, 
          here::here("data/seqone_collated_audit_data/ws136827_pdf_data.csv"),
          row.names = FALSE)

# Load INC9096 data -----------------------------------------------------------------

# This is the data for the samples that had abnormally low confidence
# in genomic instability scores due to a pipeline error (see incident INC9096).

ws140359_pdfs <- list.files(here::here("data/WS140359 data/"),
                            full.names = TRUE,
                            pattern = ".pdf")

ws140359_pdf_data <- ws140359_pdfs |> 
  map(\(ws140359_pdfs) read_seqone_report(file = ws140359_pdfs,
                                          version = "1.2")) |> 
  list_rbind() |> 
  mutate(worksheet = "WS140359_1") |> 
  filter(date == "2024-04-01")

write.csv(ws140359_pdf_data, 
          here::here("data/seqone_collated_audit_data/ws140359_pdf_data.csv"),
          row.names = FALSE)

# Collate live service PDF data -----------------------------------------------------

live_ws_excel <- read_excel(path = here::here("data/live_service_worksheets.xlsx"))

live_ws <- list(live_ws_excel$worksheet)

live_ws_filepaths_pdf <- live_ws |> 
  map(\(live_ws) find_hrd_files(live_ws, filetype = ".pdf")) |> 
  flatten()

file.copy(from = live_ws_filepaths_pdf,
          to = here::here("data/seqone_live_service_pdf_files/"))

local_live_ws_filepaths_pdf <- list.files(path = here::here("data/seqone_live_service_pdf_files/"),
                                         pattern = "*.pdf",
                                         full.names = TRUE)

collated_live_pdf_data <- local_live_ws_filepaths_pdf |>
  map(\(local_live_ws_filepaths_pdf) read_seqone_report(
    file = local_live_ws_filepaths_pdf,
    version = "1.2"
  )) |>
  list_rbind() |> 
  mutate(worksheet = ifelse(worksheet == "WS140359", "WS140359_2", worksheet))

write.csv(collated_live_pdf_data, 
          here::here("data/seqone_collated_audit_data/collated_live_pdf_data.csv"),
          row.names = FALSE)

# Collate live service CSV data -----------------------------------------------------

live_ws_filepaths_csv <- live_ws |> 
  map(\(live_ws) find_hrd_files(live_ws, filetype = ".csv")) |> 
  flatten()

file.copy(from = live_ws_filepaths_csv,
          to = here::here("data/seqone_live_service_csv_files/"))

local_live_ws_filepaths_csv <- list.files(path = here::here("data/seqone_live_service_csv_files/"),
                                          pattern = "*.csv",
                                          full.names = TRUE)

collated_live_csv_data <- local_live_ws_filepaths_csv |> 
  map(\(local_live_ws_filepaths_csv) read_seqone_csv(file = local_live_ws_filepaths_csv)) |> 
  list_rbind() |> 
  mutate(worksheet = str_extract(string = sample,
                                 pattern = "(WS[0-9]{6})_\\d{8}",
                                 group = 1))

write.csv(collated_live_csv_data, 
          here::here("data/seqone_collated_audit_data/collated_live_csv_data.csv"),
          row.names = FALSE)
