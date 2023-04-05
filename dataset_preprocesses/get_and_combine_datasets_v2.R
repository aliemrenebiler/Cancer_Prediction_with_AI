library(GEOquery)
library(Biobase)
library("biomaRt")

# Set working directory
setwd("~/Desktop")

# Set GEO data set names
gset_names_list <- list(
  "GSE13067",
  "GSE13294",
  "GSE18088",
  "GSE26682",
  "GSE35896",
  "GSE39084",
  "GSE39582",
  "GSE75316",
  "GSE92921"
)

# ------------------------------------------------------------------------------
# GET ALL DATA
# ------------------------------------------------------------------------------

# Set empty lists
gset_all_gene_names <- list()
gset_all_sample_names <- list()
gset_all_exprs <- list()
gset_all_mdata <- list()

# Set data frame for GEO set informations
gsets_info <- data.frame(matrix(nrow = 0, ncol = 3)) 
colnames(gsets_info) <- list("Dataset Name", "Sample Amount","Gene Amount")

# For every GEO set
for(gset_name in gset_names_list) {
  
  # Get the raw GEO set
  gset_raw <- getGEO(gset_name, GSEMatrix =TRUE, getGPL=FALSE)
  
  # Get the samples of platform which is named as GPL570 (if exists)
  if (length(gset_raw) > 1) {
    idx <- grep("GPL570", attr(gset_raw, "names"))
  } else {
    idx <- 1
  }
  
  # Set the necessary variables
  gset <- gset_raw[[idx]]
  gset_exprs <- as.data.frame(exprs(gset))
  gset_mdata <- pData(gset_raw[[1]])
  gset_gene_names <- rownames(gset_exprs)
  gset_sample_names <- colnames(gset_exprs)
  
  # Set GEO set information
  new_gset_info = c(
    "Dataset Name" = gset_name,
    "Sample Amount"=length(gset_sample_names),
    "Gene Amount" = length(gset_gene_names)
  )
  gsets_info[nrow(gsets_info)+1,] <- new_gset_info
    
  # Add every matrix/list to their own list
  gset_all_exprs <- append(gset_all_exprs, list(gset_exprs))
  gset_all_mdata <- append(gset_all_mdata, list(gset_mdata))
  gset_all_gene_names <- append(gset_all_gene_names, list(gset_gene_names))
  gset_all_sample_names <- append(gset_all_sample_names, list(gset_sample_names))
}

# Print all GEO set informations
message(paste0(capture.output(gsets_info), collapse = "\n"))

# ------------------------------------------------------------------------------
# SET GENE AND SAMPLE NAMES
# ------------------------------------------------------------------------------

# Intersect gene names to get mutual genes
gset_gene_names <- Reduce(intersect, gset_all_gene_names)

# Get external gene names for human (homo saphiens)
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)
hsapiens_gene_names <- getBM(
  mart=mart,
  attributes=c(
    "affy_hg_u133_plus_2",
    "external_gene_name"
  ),
  filter = "affy_hg_u133_plus_2",
  values = gset_gene_names, uniqueRows=TRUE
)

# Create a dictionary for matching the external gene names
external_to_affy_gene_names <- list()
affy_to_external_gene_names <- list()
for (i in 1:nrow(hsapiens_gene_names)) {
  # Check external gene name
  if (nchar(hsapiens_gene_names[i,"external_gene_name"]) > 0) {
    # Set gene names
    external_gene_name <- hsapiens_gene_names[i,"external_gene_name"]
    affy_gene_name <- hsapiens_gene_names[i,"affy_hg_u133_plus_2"]
    
    # Add to external to affymetrix list
    if (is.null(external_to_affy_gene_names[[external_gene_name]])){
      external_to_affy_gene_names[[external_gene_name]] <- list()
    }
    external_to_affy_gene_names[[external_gene_name]] <- append(
      external_to_affy_gene_names[[external_gene_name]],
      hsapiens_gene_names[i,"affy_hg_u133_plus_2"]
    )
    
    # Add to affymetrix to external list
    affy_to_external_gene_names[[affy_gene_name]] <- hsapiens_gene_names[i,"external_gene_name"]
  }
}

# Set gene names list as external gene names
gset_gene_names <- names(external_to_affy_gene_names)

# Union sample names to get rid of same-named ones
gset_sample_names <- Reduce(union, gset_all_sample_names)

# Set variables for output
total_sample <- length(gset_sample_names)
total_gene <- length(gset_gene_names)

message("Total Samples: ", total_sample)
message("Total Genes: ", total_gene)

# Get current date and time for file naming
current_date_time <- format(as.POSIXct(Sys.time()), format = "%y-%m-%d_%H-%M-%S")

# ------------------------------------------------------------------------------
# SET METADATA AND SAVE AS CSV
# ------------------------------------------------------------------------------

# Set combined GEO sets
gset_combined_mdata <- data.frame(matrix(nrow = length(gset_sample_names), ncol = 1)) 
colnames(gset_combined_mdata) <- list("Status")
rownames(gset_combined_mdata) <- gset_sample_names

# Set metadata values
sample_count <- 0
for (gset_mdata in gset_all_mdata) {
  
  mss_colnames <- list(
    "mmr.status:ch1",
    "msi_status:ch1",
    "microsatellite.status:ch1",
    "microsatellite status:ch1",
    "microsatellite instability (msi) status:ch1",
    "group:ch1",
    "characteristics_ch1"
  )
  
  if (length(intersect(mss_colnames, colnames(gset_mdata))) > 0) {
    # Set column name for MSS
    mss_colname <- intersect(mss_colnames, colnames(gset_mdata))[[1]]
    
    # Get sample names
    sample_names <- rownames(gset_mdata)
    
    # Set metadata values for every sample
    for (sample_name in sample_names) {
      sample_count <- sample_count + 1
      message("Getting the metadata values of sample #", sample_count, " out of ", total_sample, ": ", sample_name)
      if (grepl("MSI|dMMR", gset_mdata[sample_name, mss_colname])) {
        gset_combined_mdata[sample_name, 1] <- "MSI"
      } else {
        gset_combined_mdata[sample_name, 1] <- "MSS"
      }
    }
  }
}

message("Got all metadata values.")

# Save as csv file
write.csv(gset_combined_mdata, file=paste(current_date_time, "metadata.csv", sep = "_"))

# ------------------------------------------------------------------------------
# SET EXPRESSIONS AND SAVE AS CSV
# ------------------------------------------------------------------------------

# Set combined GEO sets
gset_combined_exprs <- data.frame(matrix(nrow = length(gset_sample_names), ncol = length(gset_gene_names))) 
colnames(gset_combined_exprs) <- gset_gene_names
rownames(gset_combined_exprs) <- gset_sample_names

# Set expression values
sample_count <- 0
for (gset_exprs in gset_all_exprs) {
  # Get sample names
  sample_names <- colnames(gset_exprs)
  
  for (sample_name in sample_names) {
    sample_count <- sample_count + 1
    message("Getting the expression values of sample #", sample_count, " out of ", total_sample, ": ", sample_name)
    for (external_gene_name in gset_gene_names) {
      
      # Get affymetrix gene names
      affy_gene_names <- external_to_affy_gene_names[[external_gene_name]]
      
      # Calculate average gene value
      gene_value <- 0
      for (affy_gene_name in affy_gene_names) {
        gene_value <- gene_value + gset_exprs[affy_gene_name, sample_name]
      }
      gene_value <- gene_value / length(affy_gene_names)
      
      # Add the average gene value
      gset_combined_exprs[sample_name, external_gene_name] <- gene_value
    }
  }
}

message("Got all expression values.")

# Save as csv file
write.csv(gset_combined_exprs, file=paste(current_date_time, "exprs.csv", sep = "_"))

# ------------------------------------------------------------------------------

message("Completed.")
