#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop(paste("4 arguments must be supplied"), call. = FALSE)
}

# Parameters
path_Rlibrary <- args[1]  # R library path
seurat_file <- args[2]    # Path to Seurat RDS file
output_dir <- args[3]     # Output directory
prefix <- args[4]         # Prefix for output files

# Libraries
.libPaths(path_Rlibrary, FALSE)
library(Seurat)

# Load Seurat object
seurat_obj <- readRDS(seurat_file)

# Filter cells
valid_cells <- !seurat_obj@meta.data$Sample_Tag %in% c("Undetermined", "Multiplet")
filtered_obj <- subset(seurat_obj, cells = rownames(seurat_obj@meta.data)[valid_cells])

# Extract subject information
singleCellSubjects <- as.character(filtered_obj@meta.data$Sample_Name)

# Save as RDA
rda_filename <- file.path(output_dir, paste0(prefix, "_single_cell_subjects.rda"))
save(singleCellSubjects, file=rda_filename)

# Save cell-subject mapping as CSV
cell_barcodes <- colnames(filtered_obj)
subjects_df <- data.frame(Cell=cell_barcodes, Subject=singleCellSubjects)
csv_filename <- file.path(output_dir, paste0(prefix, "_single_cell_subjects.csv"))
write.csv(subjects_df, file=csv_filename, row.names=FALSE)

print(paste("Single cell subjects saved to:", rda_filename, "and", csv_filename))