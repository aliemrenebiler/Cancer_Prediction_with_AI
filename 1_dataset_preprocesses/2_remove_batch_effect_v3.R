library(ggfortify)
library(limma)
library(svglite)

# ------------------------------------------------------------------------------
# CONFIGURATIONS
# ------------------------------------------------------------------------------

# Set working directory
setwd("~/Desktop")

# Set input file names
mdata_file_name <- "1_merged_mdata.csv"
exprs_file_name <- "1_merged_exprs.csv"

# Set output file names
corrected_exprs_file_name <- "2_corrected_exprs.csv"
pca_before_file_name <- "2_pca_before.svg"
pca_after_file_name <- "2_pca_after.svg"

# ------------------------------------------------------------------------------
# GET TABLES FROM CSV FILE
# ------------------------------------------------------------------------------

message("Reading CSV files...")

# Read in sample expression data
exprs_table <- read.csv(exprs_file_name, row.names = 1)

# Read in sample metadata
mdata_table <- read.csv(mdata_file_name, row.names = 1)

# Making sure the sample names in expressions matches
# to sample names in metadata in the same order
if(
  all(rownames(exprs_table) %in% rownames(mdata_table))
  && all(rownames(exprs_table) == rownames(mdata_table))
) {
  message("Sample names are matching.")
} else {
  stop("Sample names are not matching.")
}

# ------------------------------------------------------------------------------
# REMOVE BATCH EFFECT AND NORMALIZE
# ------------------------------------------------------------------------------

message('Removing batch effect...')

# Reverse column and row names, create matrix 
exprs_table_matrix <- as.matrix.data.frame(t(exprs_table))
mode(exprs_table_matrix)='numeric'

# Remove batch effect
exprs_table_corrected <- removeBatchEffect(exprs_table_matrix, mdata_table[,"Dataset"])

# Turn matrix into data frame
exprs_table_corrected <- as.data.frame(t(exprs_table_corrected))
colnames(exprs_table_corrected) <- colnames(exprs_table)
rownames(exprs_table_corrected) <- rownames(exprs_table)

# ------------------------------------------------------------------------------
# CREATE BEFORE AND AFTER PLOTS
# ------------------------------------------------------------------------------

message('Creating "before" plot...')

# Create before PCA
pca_before <- prcomp(exprs_table, scale. = TRUE)
exprs_table["Dataset"] <- mdata_table[,"Dataset"]
autoplot(pca_before, data=exprs_table, colour = 'Dataset')
ggsave(pca_before_file_name)
exprs_table <- exprs_table[ !colnames(exprs_table) %in% c('Dataset') ]

message('Creating "after" plot...')

# Create after PCA
pca_after <- prcomp(exprs_table_corrected, scale. = TRUE)
exprs_table_corrected["Dataset"] <- mdata_table[,"Dataset"]
autoplot(pca_after, data=exprs_table_corrected, colour = 'Dataset')
ggsave(pca_after_file_name)
exprs_table_corrected <- exprs_table_corrected[ !colnames(exprs_table_corrected) %in% c('Dataset')]

# Close the plots
if (length(dev.list()) > 0) {
  dev.off()
}

# ------------------------------------------------------------------------------
# SAVE AS CSV
# ------------------------------------------------------------------------------

message("Saving as CSV...")
write.csv(exprs_table_corrected, file = corrected_exprs_file_name)

# ------------------------------------------------------------------------------

message("Completed.")
