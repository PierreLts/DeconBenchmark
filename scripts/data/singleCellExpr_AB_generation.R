#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 6) {
  stop(paste("6 arguments must be supplied instead of", length(args)), call. = FALSE)
}

####### Parameter of script (ORDER IS IMPORTANT)
path_Rlibrary <- args[1] #IMPORTANT
input_data <- args[2]
output_dir <- args[3]
mapping_file <- args[4]
prefix <- args[5]
sample_filter <- args[6]  # "A", "B", or "A,B"

# Libraries
.libPaths(path_Rlibrary, FALSE) #IMPORTANT
library(Seurat)
print("CHECK: Libraries")

# Load the single-cell dataset
sc_data <- readRDS(input_data)

# Filter out "Undetermined" and "Multiplet" cells
cell_filter <- !sc_data@meta.data$Sample_Tag %in% c("Multiplet")
filtered_sc <- subset(sc_data, cells = rownames(sc_data@meta.data)[cell_filter])

# Apply A/B filtering based on sample_filter parameter
filter_type <- tolower(sample_filter)
print(paste("Applying filter:", filter_type))

# Extract sample names from Sample_Name column
sample_names <- filtered_sc@meta.data$Sample_Name

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
  stop(paste("Invalid filter type:", filter_type, "- Must be 'A', 'B', or 'A,B'"))
}

# Apply the filter
ab_filtered_sc <- subset(filtered_sc, cells = rownames(filtered_sc@meta.data)[ab_filter])

# Get count of cells for each group
a_count <- sum(grepl("\\d+A$", ab_filtered_sc@meta.data$Sample_Name))
b_count <- sum(grepl("\\d+B$", ab_filtered_sc@meta.data$Sample_Name))
print(paste("Selected cells count - A:", a_count, "B:", b_count, "Total:", ncol(ab_filtered_sc)))

# Remove unused factor levels
ab_filtered_sc@meta.data$Sample_Tag <- droplevels(ab_filtered_sc@meta.data$Sample_Tag)

# Load gene mapping file
mapping_df <- read.table(mapping_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
colnames(mapping_df) <- c("Ensembl_ID", "Gene_Name")

# Remove duplicates from mapping
mapping_df <- mapping_df[!duplicated(mapping_df$Gene_Name), ]

# Map gene names to Ensembl IDs
gene_names <- rownames(ab_filtered_sc@assays$RNA@counts)
gene_map <- setNames(mapping_df$Ensembl_ID, mapping_df$Gene_Name)
new_gene_names <- gene_map[gene_names]
new_gene_names[is.na(new_gene_names)] <- gene_names[is.na(new_gene_names)]

# Assign Ensembl IDs as rownames
rownames(ab_filtered_sc@assays$RNA@counts) <- unname(new_gene_names)

# Convert sparse matrix to dense matrix
singleCellExpr <- as.matrix(ab_filtered_sc@assays$RNA@counts)

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

# Apply duplicate handling to both datasets before finding common genes
message("Checking for duplicate genes in singleCellExpr data...")
singleCellExpr <- handle_duplicate_genes(singleCellExpr)

# Checking for duplicate columns
print('Checking for column duplicate...')
if(any(duplicated(colnames(singleCellExpr)))) {
  message("Warning: Found duplicate sample names in single-cell matrix")
  print(table(colnames(singleCellExpr))[table(colnames(singleCellExpr)) > 1])
}

# Checks
print("Preview of expression matrix (6x6):")
print(singleCellExpr[1:min(6, nrow(singleCellExpr)), 1:min(6, ncol(singleCellExpr))])
print(paste("Matrix dimensions:", paste(dim(singleCellExpr), collapse=" x ")))
print(paste("Matrix class:", class(singleCellExpr)))

# Save a small subset (50x50) of the dense matrix as CSV
max_rows <- min(50, nrow(singleCellExpr))
max_cols <- min(50, ncol(singleCellExpr))
small_dense_matrix <- singleCellExpr[1:max_rows, 1:max_cols, drop=FALSE]
small_dense_csv <- file.path(output_dir, paste0(prefix, "_singleCellExpr_", filter_label, "_50x50.csv"))
write.csv(small_dense_matrix, file = small_dense_csv)

# Save as RDA with filter label
rda_filename <- file.path(output_dir, paste0(prefix, "_singleCellExpr_", filter_label, ".rda"))
save(singleCellExpr, file = rda_filename)

# # Save full CSV
# csv_filename <- file.path(output_dir, paste0(prefix, "_singleCellExpr_", filter_label, ".csv"))
# write.csv(singleCellExpr, file = csv_filename)

print(paste("Filtered expression matrix saved to:", rda_filename))
print(paste("Filtered CSV saved to:", csv_filename))
print(paste("Small dense subset (", max_rows, "x", max_cols, ") saved to:", small_dense_csv))