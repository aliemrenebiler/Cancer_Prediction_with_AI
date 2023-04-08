library(ggfortify)
library(limma)

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
# REMOVE BATCH EFFECT
# ------------------------------------------------------------------------------

# PCA - Before
pca_before <- prcomp(exprs_table, scale. = TRUE)
exprs_table["Dataset"] <- mdata_table[,"Dataset"]
svg("pca_before.svg") # TODO: svg çalışmıyor
autoplot(pca_before, data=exprs_table, colour = 'Dataset')
exprs_table <- exprs_table[ !colnames(exprs_table) %in% c('Dataset')]

# Close the plots
dev.off()

# Reverse column and row names, create matrix 
exprs_table_matrix <- as.matrix.data.frame(t(exprs_table))
mode(exprs_table_matrix)='numeric'

# Remove batch effect
exprs_table_corrected <- removeBatchEffect(exprs_table_matrix, mdata_table[,"Dataset"])

# Turn matrix into data frame
exprs_table_corrected <- as.data.frame(t(exprs_table_corrected))
colnames(exprs_table_corrected) <- colnames(exprs_table)
rownames(exprs_table_corrected) <- rownames(exprs_table)

# PCA - After
pca_after <- prcomp(exprs_table_corrected, scale. = TRUE)
exprs_table_corrected["Dataset"] <- mdata_table[,"Dataset"]
svg("pca_after.svg") # TODO: svg çalışmıyor
autoplot(pca_after, data=exprs_table_corrected, colour = 'Dataset')
exprs_table_corrected <- exprs_table_corrected[ !colnames(exprs_table_corrected) %in% c('Dataset')]

# Close the plots
dev.off()

# TODO: Plotları yan yana bastıramıyorum

# TODO: Bunun tam amacını anlayamadım
par(mar = rep(2, 4))
png("boxplot-before.png", width = 10000,height = 5000, res=600)
boxplot(t(exprs_table), main="Original")
dev.off()
png("boxplot-after.png", width = 10000,height = 5000, res=600)
boxplot(t(exprs_table_corrected), main="Batch Corrected")
dev.off()

write.csv(exprs_table_corrected, file = corrected_exprs_file_name)
