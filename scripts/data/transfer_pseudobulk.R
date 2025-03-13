#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
  stop("5 arguments must be supplied: R_LIBRARY_PATH SOURCE_FILE MAPPING_FILE OUTPUT_DIR PREFIX SAMPLE_FILTER", call.=FALSE)
}

# Parameters
path_Rlibrary <- args[1]    # R library path
source_file <- args[2]      # Source pseudobulk counts file
mapping_file <- args[3]     # Gene name to Ensembl ID mapping file
output_dir <- args[4]       # Output directory
prefix <- args[5]           # Dataset prefix
sample_filter <- args[6]    # Sample filter: A, B, or AB

# Set library path
.libPaths(path_Rlibrary, FALSE)

# Load libraries
library(dplyr)
print("CHECK: Libraries loaded")

# Read the pseudobulk counts file
print(paste("Reading source file:", source_file))
pseudobulk_data <- read.csv(source_file, row.names=1, check.names=FALSE)
print(paste("Loaded pseudobulk data with dimensions:", paste(dim(pseudobulk_data), collapse=" x ")))

# Load gene mapping file
print(paste("Reading mapping file:", mapping_file))
mapping_df <- read.table(mapping_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
colnames(mapping_df) <- c("Ensembl_ID", "Gene_Name")

# Remove duplicates from mapping
mapping_df <- mapping_df[!duplicated(mapping_df$Gene_Name), ]
print(paste("Mapping file contains", nrow(mapping_df), "unique gene mappings"))

# Map gene names to Ensembl IDs
gene_names <- rownames(pseudobulk_data)
gene_map <- setNames(mapping_df$Ensembl_ID, mapping_df$Gene_Name)
new_gene_names <- gene_map[gene_names]
new_gene_names[is.na(new_gene_names)] <- gene_names[is.na(new_gene_names)]

# Assign Ensembl IDs as rownames
rownames(pseudobulk_data) <- unname(new_gene_names)

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

# Apply duplicate handling
message("Checking for duplicate genes in pseudobulk data...")
pseudobulk_data <- handle_duplicate_genes(as.matrix(pseudobulk_data))

# Create output directory if it doesn't exist
dest_dir <- file.path(output_dir, prefix)
dir.create(dest_dir, recursive=TRUE, showWarnings=FALSE)

# Save as RDA
rda_filename <- file.path(dest_dir, paste0(prefix, "_pseudobulk_", sample_filter, ".rda"))
pseudobulk <- pseudobulk_data  # Variable name for RDA
save(pseudobulk, file=rda_filename)

# Save as CSV
csv_filename <- file.path(dest_dir, paste0(prefix, "_pseudobulk_", sample_filter, ".csv"))
write.csv(pseudobulk_data, file=csv_filename, row.names=TRUE)

print(paste("Pseudobulk data saved to:", rda_filename, "and", csv_filename))