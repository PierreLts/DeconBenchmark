#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
  stop("5 arguments must be supplied: R_LIBRARY_PATH SOURCE_FILE MAPPING_FILE OUTPUT_DIR PREFIX", call.=FALSE)
}

# Parameters
path_Rlibrary <- args[1]    # R library path
source_file <- args[2]      # Source pseudobulk counts file
mapping_file <- args[3]     # Gene name to Ensembl ID mapping file
output_dir <- args[4]       # Output directory
prefix <- args[5]           # Dataset prefix

# Set library path
.libPaths(path_Rlibrary, FALSE)

# Load libraries
library(dplyr)
print("CHECK: Libraries loaded")

# Read the pseudobulk counts file
print(paste("Reading source file:", source_file))
pseudobulk_data <- read.csv(source_file, row.names=1, check.names=FALSE)
print(paste("Loaded pseudobulk data with dimensions:", paste(dim(pseudobulk_data), collapse=" x ")))

# Convert to matrix to ensure numeric values
pseudobulk_matrix <- as.matrix(pseudobulk_data)
print("Converted data to matrix format")

# Load gene mapping file
print(paste("Reading mapping file:", mapping_file))
mapping_df <- read.table(mapping_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
colnames(mapping_df) <- c("Ensembl_ID", "Gene_Name")

# Remove duplicates from mapping to ensure one-to-one mapping where possible
mapping_df <- mapping_df[!duplicated(mapping_df$Gene_Name), ]
print(paste("Mapping file contains", nrow(mapping_df), "unique gene mappings"))

# Extract gene names from pseudobulk data
original_genes <- rownames(pseudobulk_matrix)

# Create a named vector for mapping gene names to Ensembl IDs
gene_map <- setNames(mapping_df$Ensembl_ID, mapping_df$Gene_Name)

# Map gene names to Ensembl IDs, keeping original name if no mapping exists
ensembl_ids <- gene_map[original_genes]
ensembl_ids[is.na(ensembl_ids)] <- original_genes[is.na(ensembl_ids)]

# Create a new matrix with Ensembl IDs as row names
rownames(pseudobulk_matrix) <- ensembl_ids

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

# Apply duplicate handling to the pseudobulk data
message("Checking for duplicate genes in pseudobulk data...")
pseudobulk_matrix <- handle_duplicate_genes(pseudobulk_matrix)

# Check for column duplicates
print('Checking for column duplicates...')
if(any(duplicated(colnames(pseudobulk_matrix)))) {
  message("Warning: Found duplicate sample names in pseudobulk matrix")
  print(table(colnames(pseudobulk_matrix))[table(colnames(pseudobulk_matrix)) > 1])
}

# Create output directory if it doesn't exist
dest_dir <- file.path(output_dir, prefix)
dir.create(dest_dir, recursive=TRUE, showWarnings=FALSE)

# Save as bulk matrix for RDA (using variable name 'bulk')
bulk <- pseudobulk_matrix
rda_filename <- file.path(dest_dir, paste0(prefix, "_pseudobulk.rda"))
save(bulk, file=rda_filename)

# Save as CSV with the same format as shown in the example
csv_filename <- file.path(dest_dir, paste0(prefix, "_pseudobulk.csv"))
write.csv(bulk, file=csv_filename)

print(paste("Pseudobulk data saved to:", rda_filename, "and", csv_filename))
print(paste("Final dimensions:", paste(dim(bulk), collapse=" x ")))
print("Sample column names:")
print(head(colnames(bulk)))
print("Sample row names (Ensembl IDs):")
print(head(rownames(bulk)))