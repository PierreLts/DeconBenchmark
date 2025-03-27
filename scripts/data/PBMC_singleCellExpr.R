#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
  stop(paste("5 arguments must be supplied instead of", length(args)), call. = FALSE)
}

####### Parameter of script (ORDER IS IMPORTANT)
path_Rlibrary <- args[1] #IMPORTANT
input_data <- args[2]     # Path to PBMC Seurat object
output_dir <- args[3]     # Output directory
mapping_file <- args[4]   # Gene mapping file
prefix <- args[5]         # Prefix for output files

# Libraries
.libPaths(path_Rlibrary, FALSE) #IMPORTANT
library(Seurat)
library(Matrix)
print("CHECK: Libraries loaded")

# Load the PBMC dataset
print(paste("Loading Seurat object from:", input_data))
sc_data <- readRDS(input_data)
print("CHECK: Seurat object loaded")

# Print basic information about the object
print("=== Seurat Object Information ===")
print(sc_data)
print(paste("Number of cells:", ncol(sc_data)))
print(paste("Number of features:", nrow(sc_data)))
print(paste("Default assay:", DefaultAssay(sc_data)))

# Filter cells if needed (or use all cells)
filtered_sc <- sc_data

# Use BD_cell_type column for cell labels if available
if ("BD_cell_type" %in% colnames(filtered_sc@meta.data)) {
  print("Using BD_cell_type as cell type labels")
  cell_labels <- as.character(filtered_sc@meta.data$BD_cell_type)
} else if ("celltype" %in% colnames(filtered_sc@meta.data)) {
  print("Using celltype as cell type labels")
  cell_labels <- as.character(filtered_sc@meta.data$celltype)
} else {
  # Use active idents if no specific cell type column found
  print("No specific cell type column found, using active identities")
  cell_labels <- as.character(Idents(filtered_sc))
}

# Print cell type distribution
print("Cell type distribution:")
print(table(cell_labels, useNA = "ifany"))

# Create singleCellLabels
singleCellLabels <- cell_labels
names(singleCellLabels) <- colnames(filtered_sc)

# Create sample/subject information from orig.ident
print("Creating sample/subject information...")
if ("orig.ident" %in% colnames(filtered_sc@meta.data)) {
  print("Using 'orig.ident' as sample information")
  singleCellSubjects <- as.character(filtered_sc@meta.data$orig.ident)
  names(singleCellSubjects) <- colnames(filtered_sc)
} else {
  # Default sample name if orig.ident not available
  print("No sample information found. Using 'sample1' for all cells.")
  singleCellSubjects <- rep("sample1", length(singleCellLabels))
  names(singleCellSubjects) <- colnames(filtered_sc)
}

# Extract expression matrix
print("Extracting expression matrix...")
# Try to get counts data using Seurat's GetAssayData function
counts_data <- tryCatch({
  GetAssayData(filtered_sc, slot = "counts", assay = "RNA")
}, error = function(e) {
  # If that fails, try direct slot access
  tryCatch({
    filtered_sc@assays$RNA@counts
  }, error = function(e2) {
    # If that fails too, try accessing data in newer Seurat versions
    tryCatch({
      filtered_sc@assays$RNA@layers$counts
    }, error = function(e3) {
      stop("Could not extract count data using any method. Please check the Seurat object structure.")
    })
  })
})

# Convert to dense matrix if needed
if (inherits(counts_data, "dgCMatrix") || inherits(counts_data, "dgTMatrix")) {
  print("Converting sparse matrix to dense matrix...")
  singleCellExpr <- as.matrix(counts_data)
} else {
  singleCellExpr <- counts_data
}

print(paste("Expression matrix dimensions before gene mapping:", paste(dim(singleCellExpr), collapse=" x ")))

# === MODIFIED SECTION: Gene mapping with unmapped gene removal ===
print("Loading gene mapping file and mapping gene names to Ensembl IDs...")
if (file.exists(mapping_file)) {
  # Load gene mapping file
  mapping_df <- tryCatch({
    read.table(mapping_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
  }, error = function(e) {
    message(paste("Error reading mapping file:", e$message))
    message("Continuing without gene mapping")
    return(NULL)
  })
  
  if (!is.null(mapping_df)) {
    # Set column names if needed
    if (ncol(mapping_df) >= 2) {
      colnames(mapping_df)[1:2] <- c("Ensembl_ID", "Gene_Name")
      
      # Remove duplicates from mapping (keep first occurrence of each gene name)
      mapping_df <- mapping_df[!duplicated(mapping_df$Gene_Name), ]
      print(paste("Mapping file contains", nrow(mapping_df), "unique gene mappings"))
      
      # Map gene names to Ensembl IDs
      gene_names <- rownames(singleCellExpr)
      gene_map <- setNames(mapping_df$Ensembl_ID, mapping_df$Gene_Name)
      new_gene_names <- gene_map[gene_names]
      
      # NEW: Identify unmapped genes and remove them
      unmapped_genes <- is.na(new_gene_names)
      num_unmapped <- sum(unmapped_genes)
      num_total <- length(gene_names)
      percent_unmapped <- round(num_unmapped / num_total * 100, 2)
      
      print(paste("Total genes:", num_total))
      print(paste("Unmapped genes:", num_unmapped, "(", percent_unmapped, "%)"))
      print("Removing rows with unmapped genes...")
      
      # Filter out unmapped genes
      singleCellExpr <- singleCellExpr[!unmapped_genes, , drop=FALSE]
      new_gene_names <- new_gene_names[!unmapped_genes]
      
      # Assign new Ensembl IDs as rownames
      rownames(singleCellExpr) <- new_gene_names
      
      print(paste("Expression matrix dimensions after removing unmapped genes:", 
                 paste(dim(singleCellExpr), collapse=" x ")))
    } else {
      message("Mapping file doesn't have at least 2 columns, skipping gene mapping")
    }
  }
} else {
  print(paste("Gene mapping file not found:", mapping_file))
  print("Continuing without gene mapping")
}

# === SECTION: Handle duplicate genes function ===
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

# Apply duplicate handling to the expression matrix
print("Checking for duplicate genes in expression matrix...")
singleCellExpr <- handle_duplicate_genes(singleCellExpr)

print(paste("Final expression matrix dimensions:", paste(dim(singleCellExpr), collapse=" x ")))

# Check for column duplicates
print('Checking for column duplicates...')
if(any(duplicated(colnames(singleCellExpr)))) {
  message("Warning: Found duplicate sample names in single-cell matrix")
  print(table(colnames(singleCellExpr))[table(colnames(singleCellExpr)) > 1])
}

# Setup output directory
# Use output_dir directly without adding another prefix folder
output_dir_path <- output_dir

# Create output directory if needed
if (!dir.exists(output_dir_path)) {
  dir.create(output_dir_path, recursive = TRUE)
  print(paste("Created output directory:", output_dir_path))
}

# Save outputs - cell labels
print("Saving cell labels...")
save_path <- file.path(output_dir_path, paste0(prefix, "_singleCellLabels_AB.rda"))
save(singleCellLabels, file = save_path)

labels_df <- data.frame(
  CellBarcode = names(singleCellLabels),
  CellType = singleCellLabels,
  stringsAsFactors = FALSE
)
csv_path <- file.path(output_dir_path, paste0(prefix, "_singleCellLabels_AB.csv"))
write.csv(labels_df, file = csv_path, row.names = FALSE)

# Save outputs - expression matrix
print("Saving expression matrix...")
save_path <- file.path(output_dir_path, paste0(prefix, "_singleCellExpr_AB.rda"))
save(singleCellExpr, file = save_path)

# Save a small subset as CSV
max_rows <- min(50, nrow(singleCellExpr))
max_cols <- min(50, ncol(singleCellExpr))
small_expr <- singleCellExpr[1:max_rows, 1:max_cols]
csv_path <- file.path(output_dir_path, paste0(prefix, "_singleCellExpr_AB_50x50.csv"))
write.csv(small_expr, file = csv_path)

# Save subjects information
print("Saving cell subjects data...")
save_path <- file.path(output_dir_path, paste0(prefix, "_singleCellSubjects_AB.rda"))
save(singleCellSubjects, file = save_path)

subjects_df <- data.frame(
  CellBarcode = names(singleCellSubjects),
  Subject = singleCellSubjects,
  stringsAsFactors = FALSE
)
csv_path <- file.path(output_dir_path, paste0(prefix, "_singleCellSubjects_AB.csv"))
write.csv(subjects_df, file = csv_path, row.names = FALSE)

print("PBMC data processing completed successfully!")