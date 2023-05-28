library(GEOquery)
library(survival)
library(ggsurvfit)
library(survminer)
library(matrixStats)
library(dplyr)

# Set working directory
setwd("~/Desktop")

# ------------------------------------------------------------------------------
# CONFIGURATIONS
# ------------------------------------------------------------------------------

gset_name <- "GSE39582"

# Set file names
mdata_file_name <- "1_merged_mdata.csv"
exprs_file_name <- "3_deg_exprs.csv"
up_and_down_table_file_name <- "3_up_and_down_table.csv"

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
  message("Sample names are matching.")
} else {
  stop("Sample names are not matching.")
}

# ------------------------------------------------------------------------------
# SELECT AND COMBINE NECESSARY DATA
# ------------------------------------------------------------------------------

message("Collecting necessary data...")

# Select top 5 up and down genes
selected_up_and_down_genes <- append(
  rownames(up_and_down_table[up_and_down_table[,"type"]=="up",][1:5,]),
  rownames(up_and_down_table[up_and_down_table[,"type"]=="down",][1:5,])
)
selected_up_and_down_table <- up_and_down_table[selected_up_and_down_genes,]

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
selected_mdata_table <- mdata_table[selected_sample_names,]
selected_exprs_table <- exprs_table[selected_sample_names, selected_up_and_down_genes]

# Add new columns to metadata
selected_mdata_table <- cbind(selected_mdata_table, "OS Time"=NA)
selected_mdata_table <- cbind(selected_mdata_table, "OS Event"=NA)

# Copy overall survival times and events
selected_mdata_table[selected_sample_names,"OS Time"] <- as.integer(gset_mdata[selected_sample_names, os_time_column]) * 30
selected_mdata_table[selected_sample_names,"OS Event"] <- gset_mdata[selected_sample_names, os_event_column]

# ------------------------------------------------------------------------------
# ...
# ------------------------------------------------------------------------------

# Calculate medians of expressions for every gene
selected_up_and_down_table <- cbind(selected_up_and_down_table, "MedianExpr"=NA)
for (gene_name in selected_up_and_down_genes) {
  selected_up_and_down_table[gene_name, "MedianExpr"] <- median(as.numeric(selected_exprs_table[, gene_name]))
}

# Calculate medians of expressions for every sample
selected_mdata_table <- cbind(selected_mdata_table, "MedianExpr"=NA)
for (sample_name in selected_sample_names) {
  selected_mdata_table[sample_name, "MedianExpr"] <- median(as.numeric(selected_exprs_table[sample_name,]))
}

# Set expressions as high or low due to medians
selected_exprs_median_table <- selected_exprs_table
for (sample_name in selected_sample_names) {
  for (gene_name in selected_up_and_down_genes) {
    if (selected_exprs_table[sample_name,gene_name] > selected_up_and_down_table[gene_name,"MedianExpr"]) {
      selected_exprs_median_table[sample_name,gene_name] <- "High"
    } else {
      selected_exprs_median_table[sample_name,gene_name] <- "Low"
    }
  }
}

# Survival Fit
survival_fit <- survfit(
  Surv(
    as.numeric(selected_mdata_table[,"OS Time"]),
    as.numeric(selected_mdata_table[,"OS Event"])
  ) ~ as.numeric(selected_mdata_table[,"MedianExpr"]),
)

# Plot the result
splot <- ggsurvplot(
  survival_fit,
  data = selected_exprs_median_table,
  #conf.int = TRUE,
  surv.median.line = "hv",
  tables.y.text = FALSE,
  tables.theme = theme_cleantable(),
  xlab = "Time in Days",
  pval = T,
  # pval.method = TRUE,
  pval.coord = c(0.05, 0.05),
  risk.table = TRUE,
  legend = c(0.9, 0.2), #c("bottom")
  legend.labs = c("High","Low"),
  legend.title = c("Group"),
  title=paste0(selected_up_and_down_genes),
  fontsize = 4,
  ggtheme = theme_bw(base_family = "sans", base_size = 10)
)


print(selected_up_and_down_table[1:5,])
print(selected_mdata_table[1:5,])
print(selected_exprs_median_table)

# survfit2(Surv(as.numeric(selected_mdata_table[,"OS Time"]), as.numeric(selected_mdata_table[,"OS Event"])) ~ 1, data = selected_exprs_table) %>% 
#   ggsurvfit() +
#   labs(
#     x = "Days",
#     y = "Overall survival probability"
#   ) + 
#   add_confidence_interval() +
#   add_risktable()
