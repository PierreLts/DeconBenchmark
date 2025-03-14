#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 6) {
  stop("6 arguments must be supplied: R_LIBRARY_PATH SEURAT_OBJECT_PATH MAPPING_FILE OUTPUT_DIR PREFIX SAMPLE_FILTER", call.=FALSE)
}

# Parameters
path_Rlibrary <- args[1]  # R library path
seurat_file <- args[2]    # Seurat object file path
mapping_file <- args[3]   # Gene name to Ensembl ID mapping file
output_dir <- args[4]     # Output directory
prefix <- args[5]         # Dataset prefix
sample_filter <- args[6]  # Sample filter: A, B, or AB

# Set library path
.libPaths(path_Rlibrary, FALSE)

# Load required libraries
library(Seurat)
library(dplyr)
print("CHECK: Libraries loaded")

# Load the Seurat object
print(paste("Loading Seurat object from:", seurat_file))
seu_obj <- readRDS(seurat_file)
print("CHECK: Seurat object loaded")

# Filter cells (excluding Undetermined and Multiplet)
print("Filtering cells...")
cell_filter <- !seu_obj@meta.data$Sample_Tag %in% c("Undetermined", "Multiplet")
filtered_obj <- subset(seu_obj, cells = rownames(seu_obj@meta.data)[cell_filter])
print(paste("Retained", ncol(filtered_obj), "cells after filtering"))

# Apply A/B filtering based on sample_filter parameter
filter_type <- tolower(sample_filter)
print(paste("Applying filter:", filter_type))

# Extract sample names from Sample_Name column
sample_names <- filtered_obj@meta.data$Sample_Name

# Create A/B filter based on the parameter
if (filter_type == "a") {
  ab_filter <- grepl("\\d+A$", sample_names)
  filter_label <- "A"
} else if (filter_type == "b") {
  ab_filter <- grepl("\\d+B$", sample_names)
  filter_label <- "B"
} else if (filter_type %in% c("a,b", "b,a", "ab", "both")) {
  # Keep both A and B samples
  ab_filter <- grepl("\\d+[AB]$", sample_names)
  filter_label <- "AB"
} else {
  stop(paste("Invalid filter type:", filter_type, "- Must be 'A', 'B', or 'AB'"))
}

# Apply the filter
filtered_obj <- subset(filtered_obj, cells = rownames(filtered_obj@meta.data)[ab_filter])

# Get count of cells for each group
a_count <- sum(grepl("\\d+A$", filtered_obj@meta.data$Sample_Name))
b_count <- sum(grepl("\\d+B$", filtered_obj@meta.data$Sample_Name))
print(paste("Selected cells count - A:", a_count, "B:", b_count, "Total:", ncol(filtered_obj)))

# Extract raw counts
print("Extracting raw counts...")
counts_matrix <- as.matrix(filtered_obj@assays$RNA@counts)

# Load gene mapping file
print(paste("Reading mapping file:", mapping_file))
mapping_df <- read.table(mapping_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
colnames(mapping_df) <- c("Ensembl_ID", "Gene_Name")

# Remove duplicates from mapping to ensure one-to-one mapping where possible
mapping_df <- mapping_df[!duplicated(mapping_df$Gene_Name), ]
print(paste("Mapping file contains", nrow(mapping_df), "unique gene mappings"))

# Map gene names to Ensembl IDs
print("Mapping gene names to Ensembl IDs...")
gene_names <- rownames(counts_matrix)
gene_map <- setNames(mapping_df$Ensembl_ID, mapping_df$Gene_Name)
ensembl_ids <- gene_map[gene_names]
ensembl_ids[is.na(ensembl_ids)] <- gene_names[is.na(ensembl_ids)]  # Keep original if no mapping
rownames(counts_matrix) <- ensembl_ids

# Add Sample_Name to metadata
print("Preparing metadata for pseudobulk aggregation...")
metadata <- filtered_obj@meta.data
if (!"Sample_Name" %in% colnames(metadata)) {
  if ("orig.ident" %in% colnames(metadata)) {
    # Use orig.ident as a fallback if Sample_Name is not available
    metadata$Sample_Name <- metadata$orig.ident
    print("Using 'orig.ident' as sample grouping variable")
  } else {
    stop("Neither 'Sample_Name' nor 'orig.ident' found in metadata. Cannot determine sample grouping.")
  }
}

# Ensure cell barcodes match between counts and metadata
metadata <- metadata[match(colnames(counts_matrix), rownames(metadata)), ]

# Generate pseudobulk by summing counts across cells within each sample
print("Generating pseudobulk by aggregating cells within samples...")
sample_names_unique <- sort(unique(metadata$Sample_Name))  # Sort alphabetically
print(paste("Samples will be ordered as:", paste(sample_names_unique, collapse=", ")))
pseudobulk_counts <- matrix(0, nrow=nrow(counts_matrix), ncol=length(sample_names_unique))
rownames(pseudobulk_counts) <- rownames(counts_matrix)
colnames(pseudobulk_counts) <- sample_names_unique

# Aggregate counts for each sample
for (sample_name in sample_names_unique) {
  sample_cells <- rownames(metadata)[metadata$Sample_Name == sample_name]
  if (length(sample_cells) > 0) {
    pseudobulk_counts[, sample_name] <- rowSums(counts_matrix[, sample_cells, drop=FALSE])
  }
}

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
pseudobulk_counts <- handle_duplicate_genes(pseudobulk_counts)

# Create the bulk variable for compatibility with deconvolution pipeline
bulk <- pseudobulk_counts
print(paste("Pseudobulk matrix dimensions:", paste(dim(bulk), collapse=" x ")))

# Check for sample counts
if (ncol(bulk) == 0) {
  stop("Error: No samples found after filtering")
}
if (nrow(bulk) == 0) {
  stop("Error: No genes found after mapping and filtering")
}

# Create output directory if it doesn't exist
dest_dir <- file.path(output_dir, prefix)
dir.create(dest_dir, recursive=TRUE, showWarnings=FALSE)

# Save as RDA with filter suffix (using the naming pattern from other scripts)
rda_filename <- file.path(dest_dir, paste0(prefix, "_pseudobulk.rda"))
save(bulk, file=rda_filename)

# Save as CSV with filter suffix
csv_filename <- file.path(dest_dir, paste0(prefix, "_pseudobulk.csv"))
write.csv(bulk, file=csv_filename)

print(paste("Pseudobulk data saved to:", rda_filename, "and", csv_filename))
print(paste("Final dimensions: rows (genes):", nrow(bulk), "columns (samples):", ncol(bulk)))