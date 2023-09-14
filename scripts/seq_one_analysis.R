################################################################################
# SeqOne HRD Analysis
# joseph.shaw3@nhs.net
################################################################################

rm(list=ls())

##################################################
# Packages
##################################################

library(pdftools)
library(tidyverse)

##################################################
# Read PDFs
##################################################

location <- "~/homologous_recombination_deficiency/data/seqone_reports/"

seqone_reports <- list.files(location)

collated_info <- data.frame()

for (i in seqone_reports) {
  
  report_data <- pdftools::pdf_data(pdf = paste0(location, i))
  
  # Page 1, list 6 for text
  page1_text <- report_data[[1]][[6]]
  
  # Each value has its own coordinate within the PDF
  hrd_score <- as.numeric(page1_text[[197]])
  
  status <- page1_text[[198]]
  
  lga <- page1_text[[160]]
  
  lpc <- page1_text[[163]]
  
  ccne1 <- page1_text[[166]]
  
  rad51b <- page1_text[[169]]
  
  ncc <- page1_text[[108]]
  
  coverage <- page1_text[[123]]
  
  percent_mapping <- page1_text[[135]]
  
  sample <- page1_text[[24]]
  
  tmp_output <- data.frame("sample" = sample,
                           "hrd_score" = hrd_score,
                           "status" = status,
                           "lga" = as.numeric(lga),
                           "lpc" = as.numeric(lpc),
                           "ccne1" = as.numeric(ccne1),
                           "rad51b" = as.numeric(rad51b),
                           "ncc_percent" = as.numeric(gsub(x = ncc,
                                                           pattern = "%",
                                                           replacement = "")),
                           "coverage_x" = as.numeric(gsub(x = coverage,
                                                        pattern = "X",
                                                        replacement = "")),
                           "percent_mapping" = as.numeric(gsub(x = percent_mapping,
                                                               pattern = "%",
                                                               replacement = "")))
  
  collated_info <- rbind(collated_info, tmp_output)
  
  rm(tmp_output)
  rm(report_data)
}

##################################################
# Modify collated information
##################################################

collated_edit <- collated_info %>%
  filter(sample != "performed") %>%
  mutate(specimen_number = as.numeric(sub(x = sample,
                               pattern = "WS\\d{6}_(\\d{8})\\D{0,1}",
                               replacement = "\\1"))) %>%
  filter(specimen_number != 21012359)

##################################################
# Compare to Myriad results
##################################################

source("scripts/hrd_samples.R")

seqone_vs_myriad <- collated_edit %>%
  inner_join(sample_info, by = "specimen_number")

comparison_plot <- ggplot(seqone_vs_myriad, aes(x = gis_score, y = hrd_score)) +
  geom_point(size = 2, alpha = 0.5) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(limits = c(0,100),
                     breaks = c(0, 25, 42, 50, 75, 100)) +
  ylim(0, 1) +
  labs(x = "Myriad Genome Instability Score",
       y = "SeqOne HRD Score",
       title = "Comparison of Myriad vs SeqOne HRD Testing") +
  geom_vline(xintercept = 42, linetype = "dashed") +
  geom_hline(yintercept = 0.5, linetype = "dashed")

ggsave(plot = comparison_plot, 
       filename = paste0("comparison_plot_",
                         format(Sys.time(), "%Y%m%d_%H%M%S"),
                         ".jpg"),
       path = "~/homologous_recombination_deficiency/plots/", 
       device='jpeg',
       units = "cm",
       width = 15,
       height = 15)

##################################################