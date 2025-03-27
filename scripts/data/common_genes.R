#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
  stop("5 arguments must be supplied: <R_library_path> <input_file1> <input_file2> <output_file1> <output_file2>", call.=FALSE)
}

# Parse arguments
path_Rlibrary <- args[1]  # R library path
input_file1 <- args[2]    # First input file
input_file2 <- args[3]    # Second input file
output_file1 <- args[4]   # First output file
output_file2 <- args[5]   # Second output file

# Set library path
.libPaths(path_Rlibrary, FALSE)

# Check if input files exist
if (!file.exists(input_file1)) {
  stop(paste("Input file 1 not found:", input_file1))
}
if (!file.exists(input_file2)) {
  stop(paste("Input file 2 not found:", input_file2))
}

# Create output directories if they don't exist
dir.create(dirname(output_file1), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_file2), recursive = TRUE, showWarnings = FALSE)

# Load input files
cat("Loading first input file:", input_file1, "\n")
load(input_file1)
# The loaded data should have a variable named singleCellExpr
if (!exists("singleCellExpr")) {
  stop("Variable 'singleCellExpr' not found in the first input file")
}
expr1 <- singleCellExpr
rm(singleCellExpr)  # Remove to avoid confusion

cat("Loading second input file:", input_file2, "\n")
load(input_file2)
if (!exists("singleCellExpr")) {
  stop("Variable 'singleCellExpr' not found in the second input file")
}
expr2 <- singleCellExpr
rm(singleCellExpr)  # Remove to avoid confusion

# Get dimensions before filtering
cat("Dimensions before filtering:\n")
cat("  File 1:", paste(dim(expr1), collapse=" x "), "\n")
cat("  File 2:", paste(dim(expr2), collapse=" x "), "\n")

# Find common genes
genes1 <- rownames(expr1)
genes2 <- rownames(expr2)
common_genes <- intersect(genes1, genes2)
cat("Number of common genes:", length(common_genes), "\n")

# Filter matrices to include only common genes
expr1_filtered <- expr1[common_genes, , drop=FALSE]
expr2_filtered <- expr2[common_genes, , drop=FALSE]

# Get dimensions after filtering
cat("Dimensions after filtering:\n")
cat("  File 1:", paste(dim(expr1_filtered), collapse=" x "), "\n")
cat("  File 2:", paste(dim(expr2_filtered), collapse=" x "), "\n")

# Save filtered matrices to output files
cat("Saving filtered matrices...\n")
singleCellExpr <- expr1_filtered
save(singleCellExpr, file = output_file1)
cat("Saved to:", output_file1, "\n")

# Save a 50x50 subset as CSV for easy viewing
csv_output_file1 <- gsub("\\.rda$", "_50x50.csv", output_file1)
max_rows <- min(50, nrow(expr1_filtered))
max_cols <- min(50, ncol(expr1_filtered))
small_subset1 <- expr1_filtered[1:max_rows, 1:max_cols, drop=FALSE]
write.csv(small_subset1, file = csv_output_file1)
cat("Saved 50x50 CSV subset to:", csv_output_file1, "\n")

singleCellExpr <- expr2_filtered
save(singleCellExpr, file = output_file2)
cat("Saved to:", output_file2, "\n")

# Save a 50x50 subset as CSV for easy viewing
csv_output_file2 <- gsub("\\.rda$", "_50x50.csv", output_file2)
max_rows <- min(50, nrow(expr2_filtered))
max_cols <- min(50, ncol(expr2_filtered))
small_subset2 <- expr2_filtered[1:max_rows, 1:max_cols, drop=FALSE]
write.csv(small_subset2, file = csv_output_file2)
cat("Saved 50x50 CSV subset to:", csv_output_file2, "\n")

cat("Processing completed successfully!\n")