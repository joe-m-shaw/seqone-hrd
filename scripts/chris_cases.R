################################################################################
# Finding GEL Samples for Chris
################################################################################

rm(list=ls())

source("scripts/dlms_connection.R")

sample_columns <- sqlColumns(channel = moldb_connection,
           catalog = "MolecularDB",
           schema = "dbo",
           sqtable = "Samples")

# Add sample IDs as an input csv to avoid committing identifiers.

sample_query <- paste0("SELECT * FROM MolecularDB.dbo.Samples WHERE NGISReferralNo IN (",
                       paste(chris_cases, collapse = ", "),
                       ")")

sqlQuery(channel = moldb_connection,
         query = sample_query,
         as.is = NGISReferralNo)

################################################################################