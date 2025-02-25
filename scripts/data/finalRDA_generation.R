#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
  stop(paste("3 arguments must be supplied instead of", length(args)), call. = FALSE)
}

####### Parameter of script (ORDER IS IMPORTANT)
path_Rlibrary <- args[1] #IMPORTANT
input_data <- args[2]
output_data <- args[3]

# Libraries
.libPaths(path_Rlibrary, FALSE) #IMPORTANT

# Load data files
singleCellExpr_path <- file.path(input_data, "singleCellExpr.rda")
singleCellLabels_path <- file.path(input_data, "singleCellLabels.rda")
bulk_path <- file.path(input_data, "bulk.rda")

load(singleCellExpr_path)
load(singleCellLabels_path)
load(bulk_path)

# No need to filter cell labels - they correspond to columns in singleCellExpr
# Verify that length matches
if(length(singleCellLabels) != ncol(singleCellExpr)) {
  warning("Number of labels doesn't match number of cells after filtering")
}

# Create a list to store in .rda format
Batch1 <- list(
  bulk = bulk,
  singleCellExpr = singleCellExpr,
  singleCellLabels = singleCellLabels
)

# Checks
print(paste("bulk dimension:", paste(dim(bulk), collapse=" x ")))
print(paste("singleCellExpr dimension:", paste(dim(singleCellExpr), collapse=" x ")))
print(paste("singleCellLabels length:", length(singleCellLabels)))
print("Content of Batch1:")
print(names(Batch1))
# Save
output_file <- file.path(output_data, "Batch1.rda")
save(Batch1, file = output_file)
print(paste("Final data saved to:", output_file))