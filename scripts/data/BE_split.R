#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
  stop("3 arguments must be supplied: R_LIBRARY_PATH OUTPUT_DIR PREFIX", call.=FALSE)
}

# Parameters
path_Rlibrary <- args[1]  # R library path
output_dir <- args[2]     # Output directory
prefix <- args[3]         # Prefix for file names

# Set library path
.libPaths(path_Rlibrary, FALSE)

# Load DeconBenchmark package
library(DeconBenchmark)
print("CHECK: DeconBenchmark loaded")

# Load BloodExample data
data(BloodExample)
print("CHECK: BloodExample data loaded")
print(paste("Using prefix:", prefix))

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Extract and save bulk data
bulk <- BloodExample$bulk
save(bulk, file = file.path(output_dir, paste0(prefix, "_bulk.rda")))
print(paste("Bulk matrix saved to:", file.path(output_dir, paste0(prefix, "_bulk.rda"))))
print(paste("Bulk dimensions:", paste(dim(bulk), collapse=" x ")))

# Extract and save single cell expression data
singleCellExpr <- BloodExample$singleCellExpr
save(singleCellExpr, file = file.path(output_dir, paste0(prefix, "_singleCellExpr.rda")))
print(paste("Single cell expression matrix saved to:", file.path(output_dir, paste0(prefix, "_singleCellExpr.rda"))))
print(paste("Single cell expression dimensions:", paste(dim(singleCellExpr), collapse=" x ")))

# Extract and save single cell labels
singleCellLabels <- BloodExample$singleCellLabels
save(singleCellLabels, file = file.path(output_dir, paste0(prefix, "_singleCellLabels.rda")))
print(paste("Single cell labels saved to:", file.path(output_dir, paste0(prefix, "_singleCellLabels.rda"))))
print(paste("Number of cell labels:", length(singleCellLabels)))

# Fixed print statement
print(paste("BloodExample components have been successfully split and saved with prefix:", prefix))