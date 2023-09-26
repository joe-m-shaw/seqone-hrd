################################################################################
# Link to DNA Database on Microsoft SQL Server Studio
# joseph.shaw3@nhs.net
################################################################################

##################################################
# Packages
##################################################

library(RODBC)
library(janitor)

##################################################
# Connect to server
##################################################

# Connection setup for PC38698
moldb_connection <- RODBC::odbcConnect(dsn = "moldb")

db_tables <- sqlTables(channel = moldb_connection,
          catalog = "MolecularDB",
          schema = "dbo")

##################################################
# Function for selecting samples
##################################################

get_sample_data <- function(sample_vector) {
  
  stopifnot(length(sample_vector) >=1)
  
  sample_query <- paste0("SELECT * FROM MolecularDB.dbo.Samples WHERE LABNO IN (",
                       paste(sample_vector, collapse = ", "),
                       ")")
  
  output_data <- sqlQuery(channel = moldb_connection,
                          query = sample_query) %>%
    janitor::clean_names()
  
  
  return(output_data)
  
}

##################################################