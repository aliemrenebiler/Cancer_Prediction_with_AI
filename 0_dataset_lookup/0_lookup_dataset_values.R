library(GEOquery)

# Set working directory
setwd("~/Desktop")

# ------------------------------------------------------------------------------
# CONFIGURATIONS
# ------------------------------------------------------------------------------

# Set GEO data set name
gset_name <- "GSE39582"

# ------------------------------------------------------------------------------
# GET ALL DATA
# ------------------------------------------------------------------------------

# Get the raw GEO set
gset_raw <- getGEO(gset_name, GSEMatrix = TRUE, getGPL = FALSE)

# Get the samples of platform which is named as GPL570 (if exists)
if (length(gset_raw) > 1) {
  idx <- grep("GPL570", attr(gset_raw, "names"))
} else {
  idx <- 1
}

# Get expressions and metadata
gset <- gset_raw[[idx]]
gset_exprs <- as.data.frame(exprs(gset))
gset_mdata <- pData(gset_raw[[1]])

# ------------------------------------------------------------------------------
# PRINT ANYTHING HERE
# ------------------------------------------------------------------------------

# print(rownames(gset_mdata))
# print(length(rownames(gset_mdata)))
# print(colnames(gset_mdata))
# print(gset_mdata[1:5,])

# ------------------------------------------------------------------------------

message("Completed.")