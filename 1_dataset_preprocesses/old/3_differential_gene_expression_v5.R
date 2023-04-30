library(Biobase)
library(limma)
library(genefilter)

# Set working directory
setwd("~/Desktop")

# ------------------------------------------------------------------------------
# CONFIGURATIONS
# ------------------------------------------------------------------------------

# Set file names
mdata_file_name <- "1_merged_mdata.csv"
exprs_file_name <- "2_corrected_exprs.csv"
deg_exprs_file_name <- "3_deg_exprs.csv"
up_and_down_table_file_name <- "3_up_and_down_table.csv"

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

message("DGE analysis in progress...")

# Convert rows into columns and columns into rows as matrix
exprs_mat <- as.matrix(t(exprs_table))

# Get MSS statuses
pdata <- as.data.frame(mdata_table$Status)
rownames(pdata) <- rownames(mdata_table)
colnames(pdata) = "Label"

# Change labels as MSI or MSS
table(pdata)
pdata[pdata == "MSI-H"] <- "MSI"
pdata[pdata == "MSI-L"] <- "MSS"
table(pdata)

# Create expression set
exprs_set <- ExpressionSet(assayData = exprs_mat, phenoData = AnnotatedDataFrame(pdata))
table(pData(exprs_set)[, "Label"])

# Create the design matrix
design <- model.matrix(~0 + Label, data = pData(exprs_set))                      
head(design, 4)
colSums(design)

# Create the contrast matrix
cm <- makeContrasts(
  MSIvMSS = LabelMSI - LabelMSS,
  levels = design
)

# Fit coefficients
fit <- lmFit(exprs_set, design)

# Fit contrasts
fit2 <- contrasts.fit(fit, contrasts = cm)

# Calculate t-statistics
fit2 <- eBayes(fit2)

# Summarize results
results <- decideTests(fit2) # Default: p.value = 0.05, lfc = 0 
summary(results)

# ------------------------------------------------------------------------------
# SAVE AS CSV
# ------------------------------------------------------------------------------

message("Saving as CSV...")

# Get upregulation and downregulation
up <- which(results[,1] == 1)
down <- which(results[,1] == -1)

# Set upregulation table
up_fit <- fit2[up, ]
up_table <- topTable(up_fit[,1], number = nrow(fit2), sort.by = "logFC") 
up_table$gene <- rownames(up_table)
up_table <- up_table[order(-abs(up_table$logFC)),]
up_table$type <- "up"

# Set downregulation table
down_fit <- fit2[down, ]
down_table <- topTable(down_fit[,1], number = nrow(fit2), sort.by = "logFC") 
down_table$gene <- rownames(down_table)
down_table <- down_table[order(-abs(down_table$logFC)),]
down_table$type <- "down"

# Merge upregulation and downregulation tables and save as CSV
up_down_table <- rbind(up_table, down_table)
write.csv(up_down_table, up_and_down_table_file_name, row.names = FALSE)

# Set upregulation and downregulation gene expressions and save as CSV
exprs_dge <- exprs_table[,which(colnames(exprs_table) %in% up_down_table$gene)]
write.csv(exprs_dge, deg_exprs_file_name)

# ------------------------------------------------------------------------------

message("Completed.")













