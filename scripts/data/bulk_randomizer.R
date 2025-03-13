#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop(paste("4 arguments must be supplied instead of", length(args)), call. = FALSE)
}

####### Parameter of script (ORDER IS IMPORTANT)
path_Rlibrary <- args[1] #IMPORTANT
input_data <- args[2]
output_dir <- args[3]
prefix <- args[4]  # Prefix for output files

# Libraries
.libPaths(path_Rlibrary, FALSE) #IMPORTANT
library(Matrix)

# Load the bulk data
print(paste("Reading input data:", input_data))
feature_counts <- read.csv(input_data, row.names = 1, check.names = FALSE)

# Gene version removal
rownames(feature_counts) <- sub("\\..*", "", rownames(feature_counts))

# Create dense matrix
bulk <- as.matrix(feature_counts)

# Handle duplicate genes by summing their expression values
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
  
  # Create a new matrix for the results
  result <- matrix(0, nrow = length(unique_genes), ncol = ncol(expr_matrix))
  rownames(result) <- unique_genes
  colnames(result) <- colnames(expr_matrix)
  
  # Sum expression values for each unique gene
  for (i in seq_along(unique_genes)) {
    gene <- unique_genes[i]
    indices <- which(gene_names == gene)
    if (length(indices) > 0) {
      # Sum expression values across all instances of this gene
      result[i, ] <- colSums(expr_matrix[indices, , drop = FALSE])
    }
  }
  
  return(result)
}

# Apply duplicate handling to dataset before finding common genes
message("Checking for duplicate genes in bulk data...")
bulk <- handle_duplicate_genes(bulk)

print('Checking for column duplicates...')
if(any(duplicated(colnames(bulk)))) {
  message("Warning: Found duplicate sample names in bulk matrix")
  print(table(colnames(bulk))[table(colnames(bulk)) > 1])
}

# Now shuffle gene expression values within each sample
# This preserves the gene distribution within each sample but
# randomizes which genes have which expression values
message("Randomizing gene expression values within each sample...")
bulk_random <- bulk  # Create a copy of the original matrix

set.seed(42)  # Set seed for reproducibility
for (i in 1:ncol(bulk_random)) {
  # Shuffle expression values within this sample
  bulk_random[, i] <- sample(bulk[, i])
}


message("Randomization complete")

# Print preview of original and randomized data
cat("Original data preview:\n")
print(bulk[1:5, 1:5])

bulk <- bulk_random
cat("\nRandomized data preview:\n")
print(bulk[1:5, 1:5])

# Save as RDA
rda_filename <- file.path(output_dir, paste0(prefix, "_bulk_random.rda"))
save(bulk, file = rda_filename)

# Save as CSV
csv_filename <- file.path(output_dir, paste0(prefix, "_bulk_random.csv"))
write.csv(bulk, file = csv_filename)

print(paste("Randomized bulk matrix saved to:", rda_filename, "and", csv_filename))