#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop(paste("4 arguments must be supplied instead of", length(args)), call. = FALSE)
}

####### Parameter of script (ORDER IS IMPORTANT)
path_Rlibrary <- args[1] #IMPORTANT
input_data <- args[2]     # Input data to get dimensions from
output_dir <- args[3]     # Output directory
prefix <- args[4]         # Prefix for output files

# Libraries
.libPaths(path_Rlibrary, FALSE) #IMPORTANT
library(Matrix)

# Load original data to get dimensions
print(paste("Reading input data:", input_data))
feature_counts <- read.csv(input_data, row.names = 1, check.names = FALSE)

# Gene version removal (same as in bulk_generation.R)
rownames(feature_counts) <- sub("\\..*", "", rownames(feature_counts))

# Create a matrix with same dimensions but filled with zeros
print("Creating null matrix with same dimensions...")
bulk <- matrix(0, 
               nrow = nrow(feature_counts), 
               ncol = ncol(feature_counts))

# Keep the same row and column names
rownames(bulk) <- rownames(feature_counts)
colnames(bulk) <- colnames(feature_counts)

# Handle duplicate genes in the same way as bulk_generation.R
handle_duplicate_genes <- function(expr_matrix) {
  # Get gene names
  gene_names <- rownames(expr_matrix)
  
  # Check for duplicates
  if (!any(duplicated(gene_names))) {
    message("No duplicate genes found.")
    return(expr_matrix)
  }
  
  # Count duplicates
  dup_genes <- unique(gene_names[duplicated(gene_names)])
  dup_count <- length(dup_genes)
  message(paste("Found", dup_count, "unique genes with duplicates. Handling by summing expression values."))
  
  # Get unique gene names
  unique_genes <- unique(gene_names)
  
  # Create a new matrix for the results (still all zeros)
  result <- matrix(0, nrow = length(unique_genes), ncol = ncol(expr_matrix))
  rownames(result) <- unique_genes
  colnames(result) <- colnames(expr_matrix)
  
  return(result)
}

# Apply duplicate handling to dataset
message("Checking for duplicate genes in the null matrix...")
bulk <- handle_duplicate_genes(bulk)

# Check for column duplicates (same check as original)
print('Checking for column duplicate...')
if(any(duplicated(colnames(bulk)))) {
  message("Warning: Found duplicate sample names in bulk matrix")
  print(table(colnames(bulk))[table(colnames(bulk)) > 1])
}

# Checks
print(paste("Matrix dimensions:", paste(dim(bulk), collapse=" x ")))
print(paste("All values are zero:", all(bulk == 0)))
print(paste("Class:", class(bulk)))

# Save as RDA with _null suffix
rda_filename <- file.path(output_dir, paste0(prefix, "_bulk_null.rda"))
save(bulk, file = rda_filename)

# Save as CSV with _null suffix
csv_filename <- file.path(output_dir, paste0(prefix, "_bulk_null.csv"))
write.csv(bulk, file = csv_filename)

print(paste("Null bulk matrix saved to:", rda_filename, "and", csv_filename))