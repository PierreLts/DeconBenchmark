#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
  stop(paste("5 arguments must be supplied instead of", length(args)), call. = FALSE)
}

####### Parameter of script (ORDER IS IMPORTANT)
path_Rlibrary <- args[1] #IMPORTANT
input_data <- args[2]
output_dir <- args[3]
mapping_file <- args[4]
prefix <- args[5]  # New: prefix for output files

# Libraries
.libPaths(path_Rlibrary, FALSE) #IMPORTANT
library(Seurat)
print("CHECK: Libraries")

# Load the single-cell dataset
sc_data <- readRDS(input_data)

# Filter out "Undetermined" and "Multiplet" cells
cell_filter <- !sc_data@meta.data$Sample_Tag %in% c("Undetermined", "Multiplet")
filtered_sc <- subset(sc_data, cells = rownames(sc_data@meta.data)[cell_filter])

# Remove unused factor levels
filtered_sc@meta.data$Sample_Tag <- droplevels(filtered_sc@meta.data$Sample_Tag)

# Load gene mapping file
mapping_df <- read.table(mapping_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
colnames(mapping_df) <- c("Ensembl_ID", "Gene_Name")

# Remove duplicates from mapping
mapping_df <- mapping_df[!duplicated(mapping_df$Gene_Name), ]

# Map gene names to Ensembl IDs
gene_names <- rownames(filtered_sc@assays$RNA@counts)
gene_map <- setNames(mapping_df$Ensembl_ID, mapping_df$Gene_Name)
new_gene_names <- gene_map[gene_names]
new_gene_names[is.na(new_gene_names)] <- gene_names[is.na(new_gene_names)]

# Assign Ensembl IDs as rownames
rownames(filtered_sc@assays$RNA@counts) <- unname(new_gene_names)

# Convert sparse matrix to dense matrix
singleCellExpr <- as.matrix(filtered_sc@assays$RNA@counts)




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
print(singleCellExpr[1:6, 1:6])
print(paste("Matrix dimensions:", paste(dim(singleCellExpr), collapse=" x ")))
print(paste("Matrix class:", class(singleCellExpr)))

# Save a small subset (50x50) of the dense matrix as CSV
max_rows <- min(50, nrow(singleCellExpr))
max_cols <- min(50, ncol(singleCellExpr))
small_dense_matrix <- singleCellExpr[1:max_rows, 1:max_cols, drop=FALSE]
small_dense_csv <- file.path(output_dir, paste0(prefix, "_singleCellExpr_dense_50x50.csv"))
write.csv(small_dense_matrix, file = small_dense_csv)

# Save as RDA
rda_filename <- file.path(output_dir, paste0(prefix, "_singleCellExpr.rda"))
save(singleCellExpr, file = rda_filename)

# # Save as CSV (this might be a large file) - sparse format
# csv_filename <- file.path(output_dir, paste0(prefix, "_singleCellExpr_sparse.csv"))
# sparse_df <- as.data.frame(summary(as(singleCellExpr, "sparseMatrix")))
# colnames(sparse_df) <- c("Row", "Column", "Value")
# write.csv(sparse_df, file = csv_filename, row.names = FALSE)

# print(paste("Expression matrix saved to:", rda_filename, "and sparse representation to", csv_filename))
# print(paste("Small dense subset (", max_rows, "x", max_cols, ") saved to:", small_dense_csv))