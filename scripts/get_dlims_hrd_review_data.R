# Extract SeqOne test results from DNA Database

library(tidyverse)

DOC6590_filepath <- paste0(config::get("data_filepath"), "live_service/DOC6590/")

dbi_con <- DBI::dbConnect(
  drv = odbc::odbc(),
  dsn = "moldb")

results_tbl <- tbl(dbi_con, 
                   dbplyr::in_catalog(
                     catalog = "MolecularDB",
                     schema = "dbo",
                     table = "ResultsAccess")) |> 
  janitor::clean_names() |> 
  rename(pcrid = resultsid)

sample_tbl <- tbl(dbi_con, dbplyr::in_catalog(catalog = "MolecularDB",
                                              schema = "dbo",
                                              table = "Samples")) |> 
  janitor::clean_names()

ncc_regex <- regex(
  r"[
  (>\d{1,3}% | \d{1,3}-\d{1,3}% | \d{1,3}%-\d{1,3}% | <\d{1,3}%)
  ]",
  comments = TRUE
)

tbrca_sample_df <- sample_tbl |> 
  filter((disease == 204 | disease_2 == 204 |
           disease_3 == 204 | disease_4 == 204) &
           date_in >= "2023-12-01 00:00:00") |> 
  select(labno, disease, disease_2, disease_3, 
         disease_4, date_in, comments) |> 
  collect() |> 
  mutate(ncc_char = str_extract(string = comments, 
                           pattern = ncc_regex, 
                           group = 1))

write_csv(tbrca_sample_df, file = paste0(DOC6590_filepath,
                                        "dlims_sample_info.csv"))

tbrca_sample_vector <- tbrca_sample_df$labno

tbrca_results <- results_tbl |> 
  filter(labno %in% tbrca_sample_vector) |> 
  select(labno, pcrid, test, genotype, genotype2,
         genocomm) |> 
  collect() |> 
  mutate(worksheet = str_c("WS", pcrid)) |> 
  select(-pcrid) |> 
  relocate(worksheet, .after = labno)

seqone_results <- tbrca_results |> 
  filter(test %in% c("sWGS HRD_SeqOne", "sWGS HRD_SeqOne SSXT HS2"))

write_csv(seqone_results, file = paste0(DOC6590_filepath, 
                                        "dlims_seqone_results.csv"))

pansolid_results <- tbrca_results |> 
  filter(test %in% c(unique(grep(pattern = "seq\\span", 
                                  x = tbrca_results$test, 
                                  ignore.case = TRUE, 
                                  value = TRUE)),
                      "Proxy genotype WS for WS139270", 
                      "Qiseq re-analysis of WS139083 (CP20477)"))

write_csv(pansolid_results, file = paste0(DOC6590_filepath,
                                          "dlims_pansolid_results.csv"))
