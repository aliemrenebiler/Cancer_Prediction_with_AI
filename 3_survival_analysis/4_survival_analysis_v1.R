library(GEOquery)
library(survival)
library(survminer)
library(svglite)
library(dplyr)

# ------------------------------------------------------------------------------
# CONFIGURATIONS
# ------------------------------------------------------------------------------

# Set working directory
setwd("~/Desktop")

# Set GEO dataset name
gset_name <- "GSE39582"

# Set input file names
mdata_file_name <- "1_merged_mdata.csv"
exprs_file_name <- "3_deg_exprs.csv"
up_and_down_table_file_name <- "3_up_and_down_table.csv"

# Set output file names
survival_mdata_file_name <- "4_survival_mdata.csv"
survival_exprs_file_name <- "4_survival_exprs.csv"
survival_up_and_down_table_file_name <- "4_survival_up_and_down_table.csv"
survival_p_values_file_name <- "4_survival_analysis_p_values.csv"

# Set column names
os_time_column <- "os.delay (months):ch1"
os_event_column <- "os.event:ch1"

# ------------------------------------------------------------------------------
# GET GEO DATASET
# ------------------------------------------------------------------------------

message('Getting "', gset_name, '" dataset...')

# Get the raw GEO set
gset_raw <- getGEO(gset_name, GSEMatrix = TRUE, getGPL = FALSE)

# Get the samples of platform which is named as GPL570 (if exists)
if (length(gset_raw) > 1) {
  idx <- grep("GPL570", attr(gset_raw, "names"))
} else {
  idx <- 1
}

# Get metadata
gset <- gset_raw[[idx]]
gset_mdata <- pData(gset_raw[[1]])

# ------------------------------------------------------------------------------
# GET TABLES FROM CSV FILE
# ------------------------------------------------------------------------------

message("Reading CSV files...")

# Read in sample expression data
exprs_table <- read.csv(exprs_file_name, row.names = 1)

# Read in sample metadata
mdata_table <- read.csv(mdata_file_name, row.names = 1)

# Read in sample metadata
up_and_down_table <- read.csv(up_and_down_table_file_name)
rownames(up_and_down_table) <- up_and_down_table[,"gene"]
up_and_down_table <- up_and_down_table %>% select(-c("gene"))

# Making sure the sample names in expressions matches
# to sample names in metadata and in the same order 
if(
  all(rownames(exprs_table) %in% rownames(mdata_table))
  && all(rownames(exprs_table) == rownames(mdata_table))
  && all(rownames(up_and_down_table) %in% colnames(exprs_table))
) {
  message("All row and column names are matching.")
} else {
  stop("Some row or column names are not matching.")
}

# ------------------------------------------------------------------------------
# SELECT AND COMBINE NECESSARY DATA
# ------------------------------------------------------------------------------

message("Collecting necessary data...")

# Select top 5 up and down genes
selected_up_and_down_genes <- append(
  rownames(up_and_down_table[up_and_down_table[,"type"]=="up",][1:10,]),
  rownames(up_and_down_table[up_and_down_table[,"type"]=="down",][1:10,])
)
survival_up_and_down_table <- up_and_down_table[selected_up_and_down_genes,]

# Get selected sample names
selected_sample_names <- rownames(mdata_table[mdata_table[,"Dataset"] == gset_name,])
tmp_gset_mdata <- gset_mdata[selected_sample_names,]
selected_sample_names <- rownames(
  tmp_gset_mdata[
    tmp_gset_mdata[,os_time_column] != "N/A",
  ]
)
tmp_gset_mdata <- gset_mdata[selected_sample_names,]
selected_sample_names <- rownames(
  tmp_gset_mdata[
    tmp_gset_mdata[,os_event_column] != "N/A",
  ]
)

# Get selected parts of metadata and expressions
survival_mdata_table <- mdata_table[selected_sample_names,]
survival_exprs_table <- exprs_table[selected_sample_names, selected_up_and_down_genes]

# Add new columns to metadata
survival_mdata_table <- cbind(survival_mdata_table, "OS.Time"=NA)
survival_mdata_table <- cbind(survival_mdata_table, "OS.Event"=NA)

# Copy overall survival times and events
survival_mdata_table[selected_sample_names,"OS.Time"] <- as.integer(gset_mdata[selected_sample_names, os_time_column]) * 30
survival_mdata_table[selected_sample_names,"OS.Event"] <- gset_mdata[selected_sample_names, os_event_column]

# Calculate medians of expressions for every gene
survival_up_and_down_table <- cbind(survival_up_and_down_table, "MedianExpr"=NA)
for (gene_name in selected_up_and_down_genes) {
  survival_up_and_down_table[gene_name, "MedianExpr"] <- median(as.numeric(survival_exprs_table[, gene_name]))
}

# Set expressions as high or low due to medians
for (sample_name in selected_sample_names) {
  for (gene_name in selected_up_and_down_genes) {
    if (survival_exprs_table[sample_name,gene_name] > survival_up_and_down_table[gene_name,"MedianExpr"]) {
      survival_exprs_table[sample_name,gene_name] <- "High"
    } else {
      survival_exprs_table[sample_name,gene_name] <- "Low"
    }
  }
}

# ------------------------------------------------------------------------------
# SAVE AS CSV
# ------------------------------------------------------------------------------

message("Saving as CSV...")

write.csv(survival_mdata_table, file = survival_mdata_file_name)
write.csv(survival_exprs_table, file = survival_exprs_file_name)
write.csv(survival_up_and_down_table, file = survival_up_and_down_table_file_name)

# ------------------------------------------------------------------------------
# SURVIVAL ANALYSIS
# ------------------------------------------------------------------------------

# The method should be added in order to save survival analysis plot as SVG
grid.draw.ggsurvplot <- function(x){
  survminer:::print.ggsurvplot(x, newpage = FALSE)
}

# Create a data frame to store p-values
p_values <- data.frame(matrix(nrow = length(selected_up_and_down_genes), ncol = 2))
colnames(p_values) <- list("P.Value","Potential.Biomarker")
rownames(p_values) <- selected_up_and_down_genes

# Get time and event (censor) data
time <- as.numeric(survival_mdata_table[,"OS.Time"])
censor <- as.numeric(survival_mdata_table[,"OS.Event"])

# Run survival analysis for every gene
for(gene_name in selected_up_and_down_genes) {
  message('Survival Analysis for "', gene_name, '" gene...')
  
  gene_exprs_values <- survival_exprs_table[,gene_name]
  surv_data <- cbind.data.frame(time, censor, gene_exprs_values)
  
  # Survival Fit
  fit <- survfit(
    Surv(
      as.numeric(time),
      as.numeric(factor(censor))
    ) ~ gene_exprs_values,
    data = surv_data
  )
  
  # Save p-value and bio marker potential
  p_values[gene_name, "P.Value"] <- surv_pvalue(fit)$pval
  if (surv_pvalue(fit)$pval < 0.05) {
    p_values[gene_name, "Potential.Biomarker"] <- "Yes"
  } else {
    p_values[gene_name, "Potential.Biomarker"] <- "No"
  }
  
  # Create the plot
  survival_plot <- ggsurvplot(
    fit,
    data = surv_data,
    surv.median.line = "hv",
    tables.y.text = FALSE,
    tables.theme = theme_cleantable(),
    xlab = "Time In Days",
    pval = T,
    pval.coord = c(0.05, 0.05),
    risk.table = TRUE,
    legend = c(0.9, 0.2),
    legend.labs = c("High", "Low"),
    legend.title = c("Group"),
    title = gene_name,
    fontsize = 4,
    ggtheme = theme_bw(base_family = "sans", base_size = 10)
  )
  
  # Save the plot
  plot_save_file_name = paste0("4_gene_", gene_name, "_surv_analysis.svg")
  ggsave(plot_save_file_name, plot = survival_plot)
  
  # Close the plot(s)
  if (length(dev.list()) > 0) {
    dev.off()
  }
}

# ------------------------------------------------------------------------------
# PRINT AND SAVE RESULTS
# ------------------------------------------------------------------------------

message(paste0(capture.output(p_values), collapse = "\n"))
write.csv(p_values, file = survival_p_values_file_name)

# ------------------------------------------------------------------------------

message("Completed.")












