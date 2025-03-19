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

# Merge Seurat objects in specified order (batch 1 on top of batch 2, etc.)
cat("Merging Seurat objects...\n")
merged_seurat <- seurat_objects[[1]]
if (length(seurat_objects) > 1) {
  for (i in 2:length(seurat_objects)) {
    cat("Merging batch", i, "...\n")
    merged_seurat <- merge(
      x = merged_seurat, 
      y = seurat_objects[[i]], 
      add.cell.ids = c("", paste0("batch", i)),
      project = "MergedBatches"
    )
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