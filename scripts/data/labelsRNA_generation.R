#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
  stop(paste("3 arguments must be supplied instead of", length(args)), call. = FALSE)
}

####### Parameter of script (ORDER IS IMPORTANT)
path_Rlibrary <- args[1] #IMPORTANT
input_data <- args[2]
output_data <- args[3]

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

# Checks
print(paste("Length:", length(singleCellLabels)))
print(paste("Class of singleCellLabels (should be 'character'):", class(singleCellLabels)))

# Save only the final expression matrix
output_file <- file.path(output_data, "singleCellLabels.rda")
save(singleCellLabels, file = output_file)
print(paste("Label matrix saved to:", output_file))

