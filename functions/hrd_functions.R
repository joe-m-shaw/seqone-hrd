################################################################################
# HRD Project Functions
################################################################################

source("scripts/hrd_filepaths.R")

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

save_hrd_plot <- function(input_plot, width = 15, height = 12, dpi = 300) {
  
  # Default inputs allow for presenting a plot as half an A4 page
  
  ggsave(filename = paste0(format(Sys.time(), "%Y_%m_%d_%H_%M_%S"),
                           "_",
                           deparse(substitute(input_plot)), ".png"),
         plot = input_plot,
         device = "png",
         path = hrd_plot_path,
         units = "cm",
         width = 15,
         height = 12,
         dpi = 300)
  
}

make_individual_plot <- function(input_sample) {
  
  output_plot <- ggplot(compare_results %>%
                          filter(dlms_dna_number == input_sample), 
                        aes(x = worksheet, y = seqone_hrd_score)) +
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
               pch=21, fill=NA, size=5, colour="red", stroke=3)
  
  return(output_plot)
  
  
}

##################################################
# Table Functions
##################################################

get_sample_summary_info <- function(input_dna_no) {
  
  output <- compare_results %>%
    filter(dlms_dna_number == input_dna_no) %>%
    select(worksheet, dlms_dna_number, seqone_hrd_score,
           seqone_hrd_status, lga, lpc, ccne1, rad51b, 
           coverage.x, percent_mapping, million_reads, read_length, 
           insert_size, percent_q30, percent_aligned, percent_dups,
           myriad_gi_score, myriad_hrd_status)
  
  return(output)
  
}

##################################################