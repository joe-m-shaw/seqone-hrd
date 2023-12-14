# WS136827 Analysis

rm(list = ls())

## Packages -------------------------------------------------------------------------

library("ggpubr")
library("readxl")
library("odbc")
library("DBI")
library("dbplyr")

source("functions/hrd_functions.R")

dbi_con <- DBI::dbConnect(
  # Driver (drv) is ODBC
  drv = odbc::odbc(),
  # Data source name (dsn)
  dsn = "moldb")

# Read SeqOne PDFs ------------------------------------------------------------------

ws136827_files <- list.files("data/seq_one_reports_WS136827/",
                             full.names = TRUE)

ws136827_reports <- ws136827_files |>
  map(\(ws136827_files) read_seqone_report(
    file = ws136827_files,
    version = "1.2"
  )) |>
  list_rbind()

# Read Myriad PDFs ------------------------------------------------------------------

myriad_files <- list.files("data/myriad_reports_qiasymphony",
                           full.names = TRUE)

myriad_data <- myriad_files |>
  map(read_myriad_report) |>
  list_rbind()

# Get pathology identifiers ---------------------------------------------------------

sample_tbl <- tbl(dbi_con, dbplyr::in_catalog(catalog = "MolecularDB",
                                              schema = "dbo",
                                              table = "Samples"))

sample_numbers <- as.vector(ws136827_reports$dlms_dna_number)

sample_info <- sample_tbl |> 
  select(LABNO, PATHNO) |> 
  filter(LABNO %in% sample_numbers) |> 
  collect() |> 
  mutate(dlms_dna_number = as.numeric(LABNO))

# Join together ---------------------------------------------------------------------

extraction_details <- read_excel(path = "data/qiasymphony_validation_extraction_details.xlsx")

ws136827_joined <- extraction_details |> 
  left_join(ws136827_reports, by = "dlms_dna_number") |> 
  left_join(myriad_data, by = "nhs_number") |> 
  left_join(sample_info, by = "dlms_dna_number") |> 
  arrange(nhs_number) |> 
  mutate(check = ifelse(myriad_hrd_status == seqone_hrd_status,
                        "Status same as Myriad",
                        "Status different to Myriad"))

dual_extracted_plot <- ws136827_joined |> 
  filter(dlms_dna_number != 21003549) |> 
  ggplot(aes(x = lga, y = lpc)) +
  geom_jitter(aes(colour = seqone_hrd_status,
                  shape = extraction),
              size = 3) +
  scale_colour_manual(name = "",
                      values = c(safe_blue, safe_red, safe_grey)) +
  scale_shape_manual(name = "",
                     values = c(16, 17)) +
  geom_segment(
    data = line_df,
    mapping = aes(x = x, y = y, xend = xend, yend = yend)
  ) +
  theme_bw() +
  theme(legend.position = "bottom") +
  facet_wrap(~nhs_number) +
  labs(title = "Different extraction methods give similar SeqOne HRD results",
       subtitle = "Samples facetted by NHS number. Black line shows SeqOne positive threshold",
       x = "Large Genomic Alterations",
       y = "Loss of Parental Copy")

save_hrd_plot(dual_extracted_plot)

export_timestamp(input = ws136827_joined)

