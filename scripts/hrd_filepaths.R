################################################################################
# HRD project filepaths
################################################################################

hrd_project_path <- "S:/central shared/Genetics/Mol_Shared/Development.Team/SeqOne Homologous Recombination Deficiency Validation/HRD R script files/"

hrd_data_path <- paste0(hrd_project_path, "data/")

hrd_output_path <- paste0(hrd_project_path, "outputs/")

##################################################
# Timestamp
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