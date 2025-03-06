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

gene_names <- rownames(filtered_sc@assays$RNA@counts)
gene_map <- setNames(mapping_df$Ensembl_ID, mapping_df$Gene_Name)
new_gene_names <- gene_map[gene_names]
new_gene_names[is.na(new_gene_names)] <- gene_names[is.na(new_gene_names)]

# Detect and handle duplicate gene IDs
duplicate_genes <- new_gene_names[duplicated(new_gene_names)]
if(length(duplicate_genes) > 0) {
  cat("Detected", length(duplicate_genes), "duplicate gene IDs in single-cell data\n")
  
  # Get counts matrix
  counts_mat <- filtered_sc@assays$RNA@counts
  
  # Print all rows with duplicate gene IDs
  for(gene in unique(duplicate_genes)) {
    duplicate_indices <- which(new_gene_names == gene)
    cat("Duplicate gene:", gene, "appears", length(duplicate_indices), "times\n")
    cat("Expression values (sample of first 3 cells, max 5 values per cell):\n")
    
    for(idx in duplicate_indices) {
      # Get non-zero counts for this gene across cells (sparse matrix efficient)
      gene_counts <- counts_mat[idx, ]
      non_zero_cells <- which(gene_counts > 0)[1:min(3, sum(gene_counts > 0))]
      
      if(length(non_zero_cells) > 0) {
        cat("  Original gene:", gene_names[idx], "-> Ensembl ID:", new_gene_names[idx], "\n")
        cat("    Sample counts:", paste(sapply(non_zero_cells, function(cell) 
          paste0(colnames(counts_mat)[cell], ":", gene_counts[cell])), collapse=", "), "...\n")
      } else {
        cat("  Original gene:", gene_names[idx], "-> Ensembl ID:", new_gene_names[idx], "(All zeros)\n")
      }
    }
  }
  
  # Create a mask for first occurrences of each gene ID
  keep_mask <- !duplicated(new_gene_names)
  
  # Keep only the first occurrence of each gene
  cat("Keeping only the first occurrence of each duplicate gene ID\n")
  new_gene_names <- new_gene_names[keep_mask]
  filtered_sc@assays$RNA@counts <- filtered_sc@assays$RNA@counts[keep_mask, ]
}

# Assign Ensembl IDs as rownames
rownames(filtered_sc@assays$RNA@counts) <- unname(new_gene_names)

# Convert sparse matrix to dense matrix
singleCellExpr <- as.matrix(filtered_sc@assays$RNA@counts)

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