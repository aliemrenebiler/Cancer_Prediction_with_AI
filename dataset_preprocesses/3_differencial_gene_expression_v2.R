# Source
# https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html

library(limma)
library(edgeR)

# Set working directory
setwd("~/Desktop")

# Configs
file_date_time <- "23-04-10_19-19-18"

exprs_file_name <- paste(file_date_time, "exprs", "corrected.csv", sep = "_")
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
# PREPROCESSES
# ------------------------------------------------------------------------------

message("Preprocessing...")

# Convert rows into columns and columns into rows
exprs_table <- as.data.frame(t(exprs_table))

# Create DGE object
dge_exprs <- DGEList(exprs_table) # TODO: negative counts not allowed diyor
dge_exprs <- calcNormFactors(dge_exprs)

# Filter low-expressed genes
cutoff <- 1
drop <- which(apply(cpm(dge_exprs), 1, max) < cutoff)
dge_exprs <- dge_exprs[-drop,]

# Get sample names
sample_names <- colnames(exprs_table) 

mm <- model.matrix(mdata_table[,Status])
y <- voom(dge_exprs, mm, plot = T)

# ------------------------------------------------------------------------------
# DIFFERENTIAL GENE EXPRESSION ANALYSIS
# ------------------------------------------------------------------------------

# Fit a linear model using weighted least squares for each gene
fit <- lmFit(y, mm)

contr <- makeContrasts(groupI5.9 - groupI5.6, levels = colnames(coef(fit)))

#sfkjnkdljfnb----

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