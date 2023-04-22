library(limma)
library(genefilter)

# Set working directory
setwd("~/Desktop")

# Configs
file_date_time <- "23-04-10_19-19-18"

exprs_file_name <- paste(file_date_time, "exprs", "corrected.csv", sep = "_")
mdata_file_name <-paste(file_date_time, "metadata.csv", sep = "_")
selected_exprs_file_name <- paste(file_date_time, "exprs", "selected.csv", sep = "_")

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
# DIFFERENTIAL GENE EXPRESSION ANALYSIS
# ------------------------------------------------------------------------------

# Convert rows into columns and columns into rows
exprs_table <- as.data.frame(t(exprs_table))

# Convert the data to log2 scale
log_data <- log2(exprs_table + 1)

# Normalize the data
norm_data <- normalizeBetweenArrays(log_data, method = "quantile")

# Remove lowly expressed genes
keep <- rowSums(norm_data > 0) >= ncol(norm_data) * 0.1
norm_data <- norm_data[keep,]

# Remove highly variable genes
vars <- rowVars(norm_data)
vars[is.na(vars)] <- 0 # replace any NA values with 0
selected <- order(vars, decreasing = TRUE)[1:sum(vars > median(vars, na.rm = TRUE))]
norm_data <- norm_data[selected, ]

# Create a new column with valid R names for the Status variable
mdata_table$Status_R <- make.names(mdata_table$Status)

# Create the design matrix using the new column
design <- model.matrix(
  ~0 + Status_R,
  data = mdata_table
)

# Create the contrast matrix
contrast_matrix <- makeContrasts(
  MSI_H_vs_MSI_L = Status_RMSI.H - Status_RMSI.L,
  MSI_H_vs_MSS = Status_RMSI.H - Status_RMSS,
  MSI_L_vs_MSS = Status_RMSI.L - Status_RMSS,
  levels = colnames(design)
)

# Fit the linear model
fit <- lmFit(norm_data, design)

# Apply the eBayes function on fit
fit <- eBayes(fit)

# Perform differential expression analysis
results <- contrasts.fit(fit, contrast_matrix)

# Extract the significant genes using a cutoff of adjusted p-value < 0.05
sig_genes <- topTable(results, coef = 1, adjust.method = "BH")

# ------------------------------------------------------------------------------
# SAVE AS CSV
# ------------------------------------------------------------------------------

message("Saving as CSV...")
write.csv(sig_genes, file = selected_exprs_file_name)

# ------------------------------------------------------------------------------

message("Completed.")