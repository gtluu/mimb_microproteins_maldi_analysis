library(fuzzyjoin)
library(httr)
library(jsonlite)

# Update this with the path to your working directory containing data.
setwd('path_to_data')

# Update this with the path to the CSV file containing MALDIquant feature matrix.
# See MSV000091443.R for example MALDIquant script.
maldiquant <- read.csv('path_to_maldiquant_output_csv_file')
maldiquant <- maldiquant[, !names(maldiquant) %in% "Class"]

# Update this witht he path to the "protein.tsv" file created by FragPipe.
fragpipe <- read.table('path_to_fragpipe_output_protein.tsv_file', sep='\t', header=1)

# Get molecular weights for FragPipe protein results.
get_protein_molecular_weights <- function(uniprot_id) {
  url <- paste0("https://rest.uniprot.org/uniprotkb/", uniprot_id, ".json")
  response <- GET(url)
  if (response$status_code == 200) {
    data <- fromJSON(content(response, "text", encoding = "UTF-8"))
    return(data$sequence$molWeight)
  }
}
fragpipe$Molecular.Weight <- sapply(fragpipe$Protein.ID, get_protein_molecular_weights)
write.csv(fragpipe, 'fragpipe_protein.csv')  # Write results to prevent need for reprocessing.

# Filter to remove proteins > 30 kDa and with < 2 unique peptides.
fragpipe <- fragpipe[(fragpipe$Unique.Peptides >= 2) & (fragpipe$Molecular.Weight <= 30000),]

# Get MALDIquant masses.
maldiquant_masses <- as.numeric(sub("X", "", colnames(maldiquant)))
maldiquant_masses <- na.omit(maldiquant_masses)
maldiquant_masses <- data.frame('mz'=maldiquant_masses)

# Match MALDIquant features to FragPipe protein annotations.
fragpipe_minimal <-  fragpipe[, c('Protein.ID', 'Protein.Description')]
fragpipe_minimal$mz <- fragpipe$Molecular.Weight +  1.00782503223  # M+H weight
difference_inner_join(maldiquant_masses, fragpipe_minimal, by='mz', max_dist=0.5)
