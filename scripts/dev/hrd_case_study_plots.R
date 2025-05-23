# HRD case study data

library(tidyverse)
library(readxl)

# Read data ---------------------------------------------------------------

variants_24069447 <- read_excel(path = paste0(
  config::get("ws_filepath"),
  "WS148485/v2PANSOLID/",
  "Annotated_v2PANSOLID_WS148485_24069447_S20.xlsx"),
  sheet = "Variants_24069447") |> 
  mutate(labno = "24069447",
         case = "Case 1 (24069447)")

variants_25000823 <- read_excel(path = paste0(
  config::get("ws_filepath"),
  "WS149916/v2PANSOLID/",
  "Annotated_v2PANSOLID_WS149916_25000823_S42.xlsx"),
  sheet = "Variants_25000823") |> 
  mutate(labno = "25000823",
         case = "Case 2 (25000823)") |> 
  select(-Splice_effect)

variants_24067433 <- read_excel(path = list.files(path =  paste0(
  config::get("ws_filepath"),
  "WS148167/v2PANSOLID/"),
  full.names = TRUE,
  # Avoid adding patient name in file name
  pattern = "Annotated_v2PANSOLID_WS148167_24067433_.*.xlsx"),
  sheet = "Variants_24067433") |> 
  mutate(labno = "24067433",
         case = "Case 3 (24067433)") |> 
  select(-Splice_effect)

# Collate data ------------------------------------------------------------

variant_df <- rbind(variants_24069447, variants_25000823, variants_24067433) |> 
  janitor::clean_names() |> 
  mutate(gene = str_extract(string = variant_nomenclature,
                            pattern = "(^.*)\\sc\\..*",
                            group = 1)) |> 
  relocate(gene) |> 
  filter(!is.na(chromosome)) |> 
  filter(check_1 %in% c("TRUE", NA))

var_plot <- ggplot(variant_df, aes(x = gene, y = frequency)) +
  geom_point(aes(fill = type), shape = 21,
             size = 2, alpha = 0.7) +
  scale_fill_manual(values = c("#999999",
                               "#56B4E9",
                               "#009E73",
                               "#0072B2")) +
  geom_point(data = variant_df |> 
               filter(gene %in% c("BRCA1", "BRCA2")),
             shape = 21, size = 5, colour = "red") +
  theme_bw() +
  theme(legend.position = "bottom") +
  facet_wrap(~case, nrow = 3) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "Gene",
       y = "Variant frequency (%)",
       fill = "Variant type",
       title = "HRD case study sequence variants",
       subtitle = "Germline variants circled in red")

ggsave("hrd_case_study_plot.png", var_plot,
       units = "in", width = 6, height = 7)
