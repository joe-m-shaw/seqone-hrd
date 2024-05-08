# SeqOne LGA vs LPC Plot

library("tidyverse")
library("here")

source(here::here("functions/hrd_functions.R"))

data <- data.frame(
  lga = c(1),
  lpc = c(1)
)

new_line <- line_df |> 
  mutate(yend = ifelse(yend > 25, 25, yend))

axis_max <- 25

seqone_boundary_plot <- ggplot(data, aes(x = lga, y = lpc)) +
  theme_bw() +
  geom_segment(
    data = new_line,
    mapping = aes(x = x, y = y, xend = xend, yend = yend)) +
  labs(x = "Large Genomic Alterations", y = "Loss of Parental Copy",
       title = "SeqOne HRD Status Boundary") +
  scale_y_continuous(limits = c(0, axis_max), breaks = seq(0, axis_max, by = 1),
                     minor_breaks = FALSE) +
  scale_x_continuous(limits = c(0, axis_max), breaks = seq(0, axis_max, by = 1),
                     minor_breaks = FALSE) +
  geom_text(aes(x = 10, y = 14), vjust = 0, size = 4, label = "HRD Negative") +
  geom_text(aes(x = 21, y = 14), vjust = 0, size = 4, label = "HRD Positive")

save_hrd_plot(seqone_boundary_plot, input_width = 16, input_height = 16)
