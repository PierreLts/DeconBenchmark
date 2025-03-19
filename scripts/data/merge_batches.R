#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
  stop("At least 3 arguments must be supplied: R_LIBRARY_PATH OUTPUT_FILE BATCH1_FILE [BATCH2_FILE BATCH3_FILE ...]", call.=FALSE)
}

####### Parameter of script (ORDER IS IMPORTANT)
path_Rlibrary <- args[1]   # IMPORTANT
output_file <- args[2]     # Output path for the merged Seurat object
batch_files <- args[3:length(args)]  # Batch file paths in desired merge order

# Libraries
.libPaths(path_Rlibrary, FALSE) #IMPORTANT
library(Seurat)
library(dplyr)

cat("Will merge", length(batch_files), "Seurat objects in this order:\n")
for (i in 1:length(batch_files)) {
  cat(" ", i, ":", batch_files[i], "\n")
}

# Function to load and prepare a Seurat object with batch annotation
load_and_prepare <- function(file_path, batch_id) {
  cat("Loading batch", batch_id, "from", basename(file_path), "\n")
  seurat_obj <- readRDS(file_path)
  
  # Add batch information to metadata if not already present
  seurat_obj$batch <- paste0("batch", batch_id)
  
  # Ensure metadata columns are of correct type
  for (col in colnames(seurat_obj@meta.data)) {
    # Check if column could be a factor (character or already factor)
    if (is.character(seurat_obj@meta.data[[col]]) || is.factor(seurat_obj@meta.data[[col]])) {
      # Convert to factor to ensure consistency
      seurat_obj@meta.data[[col]] <- as.factor(seurat_obj@meta.data[[col]])
      cat("  Converted", col, "to factor\n")
    }
  }
  
  # Print basic info about the loaded object
  cat("  Features:", nrow(seurat_obj), "\n")
  cat("  Cells:", ncol(seurat_obj), "\n")
  cat("  Assays:", paste(names(seurat_obj@assays), collapse=", "), "\n")
  
  return(seurat_obj)
}

# Load all batches with batch annotation
cat("Loading Seurat objects...\n")
seurat_objects <- list()
for (i in 1:length(batch_files)) {
  seurat_objects[[i]] <- load_and_prepare(batch_files[i], i)
}

# Add batch-specific prefixes to cell names for each object
for (i in 1:length(seurat_objects)) {
  prefix <- paste0("batch", i, "_")
  cat("Adding prefix", prefix, "to cells in batch", i, "\n")
  new_cell_names <- paste0(prefix, colnames(seurat_objects[[i]]))
  seurat_objects[[i]] <- RenameCells(object = seurat_objects[[i]], new.names = new_cell_names)
}

# Merge Seurat objects
cat("Merging Seurat objects...\n")
if (length(seurat_objects) == 1) {
  merged_seurat <- seurat_objects[[1]]
} else {
  # Merge objects without additional prefixes (since we already added them)
  merged_seurat <- merge(
    x = seurat_objects[[1]],
    y = seurat_objects[-1],
    project = "MergedBatches"
  )
}

# Verify cell name consistency after merging
if (!all(rownames(merged_seurat@meta.data) == colnames(merged_seurat))) {
  warning("Cell names in metadata don't match cell names in expression matrix after merging!")
  cat("Metadata rownames mismatch count:", sum(rownames(merged_seurat@meta.data) != colnames(merged_seurat)), "\n")
}

# Make sure all metadata columns remain as factors
for (col in colnames(merged_seurat@meta.data)) {
  if (is.character(merged_seurat@meta.data[[col]]) || is.factor(merged_seurat@meta.data[[col]])) {
    merged_seurat@meta.data[[col]] <- as.factor(merged_seurat@meta.data[[col]])
    cat("Ensured", col, "is a factor in the merged object\n")
  }
}

# Basic information about the merged object
cat("\nMerged Seurat object summary:\n")
cat("Total cells:", ncol(merged_seurat), "\n")
cat("Total features:", nrow(merged_seurat), "\n")
cat("Assays:", paste(names(merged_seurat@assays), collapse=", "), "\n")
cat("Metadata columns:", paste(colnames(merged_seurat@meta.data), collapse=", "), "\n")
cat("Batch distribution:\n")
print(table(merged_seurat$batch))

# Save the merged Seurat object
cat("Saving merged Seurat object to", output_file, "\n")
saveRDS(merged_seurat, file = output_file)

cat("Merge complete!\n")