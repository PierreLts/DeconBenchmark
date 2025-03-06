#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop(paste("4 arguments must be supplied instead of", length(args)), call. = FALSE)
}

####### Parameter of script (ORDER IS IMPORTANT)
path_Rlibrary <- args[1] #IMPORTANT
input_data <- args[2]
output_dir <- args[3]
prefix <- args[4]  # New: prefix for output files

# Libraries
.libPaths(path_Rlibrary, FALSE) #IMPORTANT
library(Seurat)

# Load
seurat_obj <- readRDS(input_data)
# Filter out "Undetermined" and "Multiplet" cells
cell_filter <- !seurat_obj@meta.data$Sample_Tag %in% c("Undetermined", "Multiplet")
filtered_obj <- subset(seurat_obj, cells = rownames(seurat_obj@meta.data)[cell_filter])
# Extract cell labels
cell_labels <- filtered_obj@meta.data$Cell_Type_Experimental
# Character vector
singleCellLabels <- as.character(cell_labels)
names(singleCellLabels) <- colnames(filtered_obj)

# Checks
print(paste("Length:", length(singleCellLabels)))
print(paste("Class of singleCellLabels (should be 'character'):", class(singleCellLabels)))

# Save as RDA
rda_filename <- file.path(output_dir, paste0(prefix, "_singleCellLabels.rda"))
save(singleCellLabels, file = rda_filename)

# Save as CSV
csv_filename <- file.path(output_dir, paste0(prefix, "_singleCellLabels.csv"))
write.csv(data.frame(CellBarcode = names(singleCellLabels), CellType = singleCellLabels), 
          file = csv_filename, row.names = FALSE)

print(paste("Cell labels saved to:", rda_filename, "and", csv_filename))