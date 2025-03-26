#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
  stop(paste("5 arguments must be supplied instead of", length(args)), call. = FALSE)
}

####### Parameter of script (ORDER IS IMPORTANT)
path_Rlibrary <- args[1] #IMPORTANT
input_data <- args[2]     # Path to PBMC Seurat object
output_dir <- args[3]     # Output directory
mapping_file <- args[4]   # Gene mapping file (not used in simplified version)
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

print(paste("Expression matrix dimensions:", paste(dim(singleCellExpr), collapse=" x ")))

# Setup output directories
output_dir_base <- output_dir
output_dir_prefix <- file.path(output_dir_base, prefix)

# Create output directories if needed
if (!dir.exists(output_dir_base)) {
  dir.create(output_dir_base, recursive = TRUE)
}
if (!dir.exists(output_dir_prefix)) {
  dir.create(output_dir_prefix, recursive = TRUE)
}

# Save outputs - cell labels
print("Saving cell labels...")
save_path <- file.path(output_dir_prefix, paste0(prefix, "_singleCellLabels_AB.rda"))
save(singleCellLabels, file = save_path)

labels_df <- data.frame(
  CellBarcode = names(singleCellLabels),
  CellType = singleCellLabels,
  stringsAsFactors = FALSE
)
csv_path <- file.path(output_dir_prefix, paste0(prefix, "_singleCellLabels_AB.csv"))
write.csv(labels_df, file = csv_path, row.names = FALSE)

# Save outputs - expression matrix
print("Saving expression matrix...")
save_path <- file.path(output_dir_prefix, paste0(prefix, "_singleCellExpr_AB.rda"))
save(singleCellExpr, file = save_path)

# Save a small subset as CSV
max_rows <- min(50, nrow(singleCellExpr))
max_cols <- min(50, ncol(singleCellExpr))
small_expr <- singleCellExpr[1:max_rows, 1:max_cols]
csv_path <- file.path(output_dir_prefix, paste0(prefix, "_singleCellExpr_AB_50x50.csv"))
write.csv(small_expr, file = csv_path)

# Save subjects information
print("Saving cell subjects data...")
save_path <- file.path(output_dir_prefix, paste0(prefix, "_singleCellSubjects_AB.rda"))
save(singleCellSubjects, file = save_path)

subjects_df <- data.frame(
  CellBarcode = names(singleCellSubjects),
  Subject = singleCellSubjects,
  stringsAsFactors = FALSE
)
csv_path <- file.path(output_dir_prefix, paste0(prefix, "_singleCellSubjects_AB.csv"))
write.csv(subjects_df, file = csv_path, row.names = FALSE)

print("PBMC data processing completed successfully!")