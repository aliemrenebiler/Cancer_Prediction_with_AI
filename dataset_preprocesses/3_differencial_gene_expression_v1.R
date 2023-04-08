# Video Sources:
# https://www.youtube.com/watch?v=UFB993xufUU
# https://www.youtube.com/watch?v=OzNzO8qwwp0

library(DESeq2)
library(tidyverse)
library(airway)

# Set working directory
setwd("~/Desktop")

# Configs
file_date_time <- "23-04-05_19-39-00"

exprs_file_name <- paste(file_date_time, "exprs.csv", sep = "_")
mdata_file_name <-paste(file_date_time, "metadata.csv", sep = "_")

# ------------------------------------------------------------------------------
# PREPARING COUNT DATA
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
# CONSTRUCT A DESEQ DATASET OBJECT
# ------------------------------------------------------------------------------

# Convert rows into columns and columns into rows
exprs_table <- as.data.frame(t(exprs_table))

# Construct DESeq Object
dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(exprs_table),
  colData = mdata_table,
  design = ~ Status
)

# TODO: Is this neccessary? How we can do that?
# # Pre-Filtering: Removing rows with low gene expression,
# # keeping rows that have at least 10 reads total
# keep <- rowSums(counts(dds)) >= 10
# dds <- dds[keep,]

# Set the factor level
dds$Status <- relevel(dds$Status, ref = "MSI-H")

# ------------------------------------------------------------------------------
# RUNNING DIFFERENTIAL GENE EXPRESSION ANALYSIS AND EXPLORE RESULTS
# ------------------------------------------------------------------------------

# Run DESeq
dds <- DESeq(dds)
res <- results(dds)

# Explore Results
summary(res)

res0.01 <- results(dds, alpha = 0.01)
summary(res0.01)

# Contrasts
resultsNames(dds)

# TODO: What is this?
# # e.g.: treated_4hrs, treated_8hrs, untreated
# results(dds, contrast = c("dexamethasone", "treated_4hrs", "untreated"))

# MA plot
plotMA(res)