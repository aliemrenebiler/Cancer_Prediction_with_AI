library(Biobase)
library(limma)
library(genefilter)

# Set working directory
setwd("~/Desktop")

# Configs
file_date_time <- "23-04-10_19-19-18"

exprs_file_name <- paste(file_date_time, "exprs", "corrected.csv", sep = "_")
mdata_file_name <-paste(file_date_time, "metadata.csv", sep = "_")
exprs_after_dge_file_name <- paste(file_date_time, "exprs", "after", "dge.csv", sep = "_")

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

up <- which(results[,1] == 1)
upGenes <- fit2[up, ]
upTable <- topTable(upGenes[,1], number = nrow(fit2), sort.by = "logFC") 
upTable$gene <- rownames(upTable)
upTable <- upTable[order(-abs(upTable$logFC)),]
rn_up = upTable$gene

write.csv(rn_up, "up-genes.csv",row.names = FALSE)
write.csv(upTable, "upTable.csv",row.names = FALSE)

down <- which(results[,1] == -1)
downGenes <- fit2[down, ]
downTable <- topTable(downGenes[,1], number = nrow(fit2), sort.by = "logFC") 
downTable$gene <- rownames(downTable)
downTable <- downTable[order(-abs(downTable$logFC)),]
rn_down = downTable$gene

write.csv(rn_down, "down-genes.csv",row.names = FALSE)
write.csv(downTable, "downTable.csv",row.names = FALSE)

dge_genes <- c(rn_up,rn_down)
write.csv(dge_genes, "dgeGenes.csv",row.names = FALSE)

exprs_dge<-exprs_table[,which(colnames(exprs_table) %in% dge_genes)]
write.csv(exprs_dge, exprs_after_dge_file_name)

# ------------------------------------------------------------------------------

message("Completed.")













