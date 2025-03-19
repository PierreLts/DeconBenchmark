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
# Also save as CSV
write.csv(bulk, file = file.path(output_dir, paste0(prefix, "_bulk.csv")))
print(paste("Bulk matrix saved to:", file.path(output_dir, paste0(prefix, "_bulk.rda")), "and", file.path(output_dir, paste0(prefix, "_bulk.csv"))))
print(paste("Bulk dimensions:", paste(dim(bulk), collapse=" x ")))

# Extract and save single cell expression data
singleCellExpr <- BloodExample$singleCellExpr
save(singleCellExpr, file = file.path(output_dir, paste0(prefix, "_singleCellExpr.rda")))

# Save a small subset (50x50) as CSV - the full matrix would be too large
max_rows <- min(50, nrow(singleCellExpr))
max_cols <- min(50, ncol(singleCellExpr))
small_expr_matrix <- singleCellExpr[1:max_rows, 1:max_cols, drop=FALSE]
# Save the subset with a descriptive name
write.csv(small_expr_matrix, file = file.path(output_dir, paste0(prefix, "_singleCellExpr_50x50.csv")))

print(paste("Single cell expression matrix saved to:", file.path(output_dir, paste0(prefix, "_singleCellExpr.rda"))))
print(paste("Small subset (50x50) saved as CSV to:", file.path(output_dir, paste0(prefix, "_singleCellExpr_50x50.csv"))))
print(paste("Single cell expression dimensions:", paste(dim(singleCellExpr), collapse=" x ")))

# Extract and save single cell labels
singleCellLabels <- BloodExample$singleCellLabels
save(singleCellLabels, file = file.path(output_dir, paste0(prefix, "_singleCellLabels.rda")))
# Create a data frame for the labels and save as CSV
# Check if names exist, if not create them
if (is.null(names(singleCellLabels))) {
  # If no names, create cell identifiers
  cell_ids <- paste0("cell_", 1:length(singleCellLabels))
} else {
  cell_ids <- names(singleCellLabels)
}
labels_df <- data.frame(
  CellID = cell_ids,
  CellType = singleCellLabels,
  stringsAsFactors = FALSE
)
write.csv(labels_df, file = file.path(output_dir, paste0(prefix, "_singleCellLabels.csv")), row.names = FALSE)
print(paste("Single cell labels saved to:", file.path(output_dir, paste0(prefix, "_singleCellLabels.rda")), "and", file.path(output_dir, paste0(prefix, "_singleCellLabels.csv"))))
print(paste("Number of cell labels:", length(singleCellLabels)))

# Create a ground truth proportions file based on the cell type distribution
cell_type_counts <- table(singleCellLabels)
proportions <- prop.table(cell_type_counts)
proportions_df <- data.frame(
  CellType = names(proportions),
  Proportion = as.numeric(proportions),
  stringsAsFactors = FALSE
)
write.csv(proportions_df, file = file.path(output_dir, paste0(prefix, "_GT_proportions.csv")), row.names = FALSE)

# Also save as RDA in the format expected by benchmarking scripts
P <- matrix(as.numeric(proportions), nrow = 1, dimnames = list("true_proportions", names(proportions)))
groundTruth <- list(P = P)
save(groundTruth, file = file.path(output_dir, paste0(prefix, "_GT_proportions.rda")))

print(paste("Ground truth proportions saved to:", file.path(output_dir, paste0(prefix, "_GT_proportions.rda")), "and", file.path(output_dir, paste0(prefix, "_GT_proportions.csv"))))

# Fixed print statement
print(paste("BloodExample components have been successfully split and saved with prefix:", prefix))