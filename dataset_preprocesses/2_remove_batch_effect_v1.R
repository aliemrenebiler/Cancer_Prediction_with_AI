library(ggfortify)
library(limma)
library(svglite)

# Set working directory
setwd("~/Desktop")

# Configs
file_date_time <- "23-04-05_19-39-00"

exprs_file_name <- paste(file_date_time, "exprs.csv", sep = "_")
mdata_file_name <-paste(file_date_time, "metadata.csv", sep = "_")
corrected_exprs_file_name <- paste(file_date_time, "exprs", "corrected.csv", sep = "_")

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

message('Normalizing...')

# Normalize the data frame between 0 and 1
exprs_table_normalized = (exprs_table_corrected - min(exprs_table_corrected)) / (max(exprs_table_corrected) - min(exprs_table_corrected))

# ------------------------------------------------------------------------------
# CREATE BEFORE AND AFTER PLOTS
# ------------------------------------------------------------------------------

# TODO: Plotları yan yana bastıramıyorum

message('Creating "before" plot...')

# Create before PCA
pca_before <- prcomp(exprs_table, scale. = TRUE)
exprs_table["Dataset"] <- mdata_table[,"Dataset"]
autoplot(pca_before, data=exprs_table, colour = 'Dataset')
ggsave("pca_before.svg")
exprs_table <- exprs_table[ !colnames(exprs_table) %in% c('Dataset')]

# Close the plots
dev.off()

message('Creating "after" plot...')

# Create after PCA
pca_after <- prcomp(exprs_table_normalized, scale. = TRUE)
exprs_table_normalized["Dataset"] <- mdata_table[,"Dataset"]
autoplot(pca_after, data=exprs_table_normalized, colour = 'Dataset')
ggsave("pca_after.svg")
exprs_table_normalized <- exprs_table_normalized[ !colnames(exprs_table_normalized) %in% c('Dataset')]

# Close the plots
dev.off()

message("Creating other plots...")

# TODO: Bunun tam amacını anlayamadım
par(mar = rep(2, 4))
png("boxplot_before.png", width = 10000,height = 5000, res=600)
boxplot(t(exprs_table), main="Original")
dev.off()
png("boxplot_after.png", width = 10000,height = 5000, res=600)
boxplot(t(exprs_table_normalized), main="Batch Corrected and Normalized")
dev.off()

# ------------------------------------------------------------------------------
# SAVE AS CSV
# ------------------------------------------------------------------------------

message("Saving as CSV...")
write.csv(exprs_table_normalized, file = corrected_exprs_file_name)

# ------------------------------------------------------------------------------

message("Completed.")