#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop(paste("4 arguments must be supplied instead of", length(args)), call. = FALSE)
}

####### Parameter of script (ORDER IS IMPORTANT)
path_Rlibrary <- args[1] #IMPORTANT
input_dir <- args[2]
output_dir <- args[3]
prefix <- args[4]  # New: prefix for output files

# Libraries
.libPaths(path_Rlibrary, FALSE) #IMPORTANT

# Load data files
singleCellExpr_path <- file.path(input_dir, paste0(prefix, "_singleCellExpr.rda"))
singleCellLabels_path <- file.path(input_dir, paste0(prefix, "_singleCellLabels.rda"))
bulk_path <- file.path(input_dir, paste0(prefix, "_bulk.rda"))

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
output_file <- file.path(output_dir, paste0(prefix, "_Batch1.rda"))
save(Batch1, file = output_file)
print(paste("Final data saved to:", output_file))

# CSV file of Batch1 structure
csv_info <- data.frame(
  Component = c("bulk", "singleCellExpr", "singleCellLabels"),
  RowCount = c(nrow(bulk), nrow(singleCellExpr), 1),
  ColumnCount = c(ncol(bulk), ncol(singleCellExpr), length(singleCellLabels)),
  Description = c("Bulk RNA-seq data", "Single-cell expression data", "Cell type labels")
)
csv_filename <- file.path(output_dir, paste0(prefix, "_Batch1_structure.csv"))
write.csv(csv_info, file = csv_filename, row.names = FALSE)
print(paste("Batch1 structure info saved to:", csv_filename))