#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
  stop(paste("5 arguments must be supplied instead of", length(args)), call. = FALSE)
}

####### Parameter of script (ORDER IS IMPORTANT)
path_Rlibrary <- args[1] 
seurat_file <- args[2]  
output_dir <- args[3]  
prefix <- args[4]  
sample_filter <- args[5]  # New parameter

# Libraries
.libPaths(path_Rlibrary, FALSE) #IMPORTANT
library(Seurat)
library(reshape2)

# Load the Seurat object
sc_data <- readRDS(seurat_file)

# Filter out "Undetermined" and "Multiplet" cells
valid_cells <- !sc_data@meta.data$Sample_Tag %in% c("Undetermined", "Multiplet")
filtered_sc <- subset(sc_data, cells = rownames(sc_data@meta.data)[valid_cells])

# Get the valid sample names
valid_sample_names <- unique(filtered_sc@meta.data$Sample_Name)

# Create a list to store the ground truth proportions for each sample
sample_proportions <- list()

# Calculate the ground truth proportions per sample
for (sample_name in valid_sample_names) {
  # Extract cells for this sample
  sample_cells <- rownames(filtered_sc@meta.data)[filtered_sc@meta.data$Sample_Name == sample_name]
  
  # Get cell types for these cells
  sample_cell_types <- filtered_sc@meta.data$Cell_Type_Experimental[filtered_sc@meta.data$Sample_Name == sample_name]
  
  # Calculate proportions
  cell_counts <- table(sample_cell_types)
  proportions <- as.numeric(cell_counts) / sum(cell_counts)
  
  # Store as a named vector
  sample_proportions[[sample_name]] <- setNames(proportions, names(cell_counts))
}

# Create a matrix with samples as rows and cell types as columns
# First, get all unique cell types across all samples
all_cell_types <- unique(unlist(lapply(sample_proportions, names)))

# Initialize an empty matrix
P <- matrix(0, 
            nrow = length(sample_proportions), 
            ncol = length(all_cell_types),
            dimnames = list(names(sample_proportions), all_cell_types))

# Fill in the matrix with the proportions
for (sample_name in names(sample_proportions)) {
  sample_props <- sample_proportions[[sample_name]]
  P[sample_name, names(sample_props)] <- sample_props
}

# Create a list to store the ground truth
groundTruth <- list(P = P)

# Save as RDA
rda_filename <- file.path(output_dir, paste0(prefix, "_GT_proportions_per_sample.rda"))
save(groundTruth, file = rda_filename)

# Save as CSV
csv_filename <- file.path(output_dir, paste0(prefix, "_GT_proportions_per_sample.csv"))
# Convert matrix to long format for CSV
gt_df <- as.data.frame(groundTruth$P)
gt_df$Sample <- rownames(gt_df)
gt_long <- melt(gt_df, id.vars = "Sample", variable.name = "CellType", value.name = "Proportion")
write.csv(gt_long, file = csv_filename, row.names = FALSE)

print(paste("Per-sample ground truth saved to:", rda_filename, "and", csv_filename))