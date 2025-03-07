#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop(paste("4 arguments must be supplied instead of", length(args)), call. = FALSE)
}

####### Parameter of script (ORDER IS IMPORTANT)
path_Rlibrary <- args[1] #IMPORTANT
seurat_file <- args[2]  # Path to Seurat RDS file (not Batch1.rda)
output_dir <- args[3]  # Output directory
prefix <- args[4]  # Prefix for output files

# Libraries
.libPaths(path_Rlibrary, FALSE) #IMPORTANT
library(Seurat)

# Load the Seurat object
sc_data <- readRDS(seurat_file)

# Filter out "Undetermined" and "Multiplet" cells
valid_cells <- !sc_data@meta.data$Sample_Tag %in% c("Undetermined", "Multiplet")
filtered_sc <- subset(sc_data, cells = rownames(sc_data@meta.data)[valid_cells])

# Get all cell types and calculate overall proportions
cell_types <- filtered_sc@meta.data$Cell_Type_Experimental
cell_counts <- table(cell_types)
proportions <- as.numeric(cell_counts) / sum(cell_counts)

# Create a matrix with cell types as columns
P <- matrix(proportions, nrow = 1, dimnames = list("true_proportions", names(cell_counts)))

# Create a list to store the ground truth
groundTruth <- list(P = P)

# Print information for verification
print("Ground truth proportions:")
print(groundTruth$P)
print(paste("Total cell types:", length(unique(cell_types))))
print(paste("Sum of proportions:", sum(proportions)))

# Save as RDA
rda_filename <- file.path(output_dir, paste0(prefix, "_ground_truth_proportions.rda"))
save(groundTruth, file = rda_filename)

# Save as CSV
csv_filename <- file.path(output_dir, paste0(prefix, "_ground_truth_proportions.csv"))
gt_df <- as.data.frame(t(groundTruth$P))
gt_df$CellType <- rownames(gt_df)
colnames(gt_df) <- c("Proportion", "CellType")
gt_df <- gt_df[, c("CellType", "Proportion")]
write.csv(gt_df, file = csv_filename, row.names = FALSE)

print(paste("Ground truth saved to:", rda_filename, "and", csv_filename))