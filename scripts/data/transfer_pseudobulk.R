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

# Load gene mapping file
print(paste("Reading mapping file:", mapping_file))
mapping_df <- read.table(mapping_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
colnames(mapping_df) <- c("Ensembl_ID", "Gene_Name")

# Remove duplicates from mapping
mapping_df <- mapping_df[!duplicated(mapping_df$Gene_Name), ]
print(paste("Mapping file contains", nrow(mapping_df), "unique gene mappings"))

# Create a matrix with the original data
pseudobulk_matrix <- as.matrix(pseudobulk_data)
gene_names <- rownames(pseudobulk_matrix)

# Map gene names to Ensembl IDs without setting as row names yet
gene_map <- setNames(mapping_df$Ensembl_ID, mapping_df$Gene_Name)
new_gene_names <- gene_map[gene_names]
new_gene_names[is.na(new_gene_names)] <- gene_names[is.na(new_gene_names)]

# Create a new data frame with Ensembl IDs and original data
# Important: Don't set row names yet to avoid duplicate error
mapped_data <- data.frame(
  Ensembl_ID = new_gene_names,
  pseudobulk_matrix,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Handle duplicate genes by summing their expression values
handle_duplicate_genes <- function(data_frame) {
  # First column contains Ensembl IDs
  ensembl_ids <- data_frame[, 1]
  expr_data <- data_frame[, -1, drop = FALSE]
  
  # Check for duplicates
  if (!any(duplicated(ensembl_ids))) {
    message("No duplicate genes found.")
    result_matrix <- as.matrix(expr_data)
    rownames(result_matrix) <- ensembl_ids
    return(result_matrix)
  }
  
  # Count duplicates
  dup_genes <- unique(ensembl_ids[duplicated(ensembl_ids)])
  dup_count <- length(dup_genes)
  message(paste("Found", dup_count, "unique genes with duplicates. Handling by summing expression values."))
  
  # Get unique gene names
  unique_genes <- unique(ensembl_ids)
  
  # Create a new matrix for the results
  result <- matrix(0, nrow = length(unique_genes), ncol = ncol(expr_data))
  rownames(result) <- unique_genes
  colnames(result) <- colnames(expr_data)
  
  # Sum expression values for each unique gene
  for (i in seq_along(unique_genes)) {
    gene <- unique_genes[i]
    indices <- which(ensembl_ids == gene)
    if (length(indices) > 0) {
      # Sum expression values across all instances of this gene
      result[i, ] <- colSums(expr_data[indices, , drop = FALSE])
    }
  }
  
  return(result)
}

# Apply duplicate handling
message("Checking for duplicate genes in pseudobulk data...")
pseudobulk_matrix <- handle_duplicate_genes(mapped_data)

# Create output directory if it doesn't exist
dest_dir <- file.path(output_dir, prefix)
dir.create(dest_dir, recursive=TRUE, showWarnings=FALSE)

# Save as RDA
rda_filename <- file.path(dest_dir, paste0(prefix, "_pseudobulk.rda"))
bulk <- pseudobulk_matrix  # Variable name for RDA
save(bulk, file=rda_filename)

# Save as CSV
csv_filename <- file.path(dest_dir, paste0(prefix, "_pseudobulk.csv"))
write.csv(pseudobulk_matrix, file=csv_filename, row.names=TRUE)

print(paste("Pseudobulk data saved to:", rda_filename, "and", csv_filename))