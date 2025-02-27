#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
  stop(paste("3 arguments must be supplied instead of", length(args)), call. = FALSE)
}

####### Parameter of script (ORDER IS IMPORTANT)
path_Rlibrary <- args[1] #IMPORTANT
seurat_file <- args[2]  # Path to Seurat RDS file
output_file <- args[3]  # Full path to output ground truth file

# Libraries
print("Setting library paths")
.libPaths(path_Rlibrary, FALSE) #IMPORTANT
print(.libPaths())

print("Loading Seurat library")
library(Seurat)

# Check if file exists
print(paste("Checking if file exists:", seurat_file))
if (!file.exists(seurat_file)) {
  stop(paste("Seurat file not found:", seurat_file))
}

# Load the Seurat object
print(paste("Loading Seurat object from:", seurat_file))
sc_data <- readRDS(seurat_file)
print("Successfully loaded Seurat object")
print(paste("Number of cells:", ncol(sc_data)))
print(paste("Number of features:", nrow(sc_data)))

# Filter out "Undetermined" and "Multiplet" cells
print("Filtering out Undetermined and Multiplet cells")
valid_cells <- !sc_data@meta.data$Sample_Tag %in% c("Undetermined", "Multiplet")
filtered_sc <- subset(sc_data, cells = rownames(sc_data@meta.data)[valid_cells])
print(paste("Cells after filtering:", ncol(filtered_sc)))

# Get the valid sample names
valid_sample_names <- unique(filtered_sc@meta.data$Sample_Name)
print("Valid sample names:")
print(valid_sample_names)

# Create a list to store the ground truth proportions for each sample
sample_proportions <- list()

# Calculate the ground truth proportions per sample
print("Calculating proportions per sample...")
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
  
  # Print information for verification
  print(paste("Sample:", sample_name))
  print(paste("  Total cells:", length(sample_cells)))
  print(paste("  Cell types found:", length(cell_counts)))
  print(paste("  Sum of proportions:", sum(proportions)))
}

# Create a matrix with samples as rows and cell types as columns
print("Creating ground truth proportion matrix...")
# First, get all unique cell types across all samples
all_cell_types <- unique(unlist(lapply(sample_proportions, names)))
print(paste("Total unique cell types across all samples:", length(all_cell_types)))

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

# Print the final ground truth matrix summary
print("Ground truth proportions matrix dimensions:")
print(dim(groundTruth$P))
print("First few rows and columns:")
print(head(groundTruth$P[, 1:min(5, ncol(groundTruth$P))]))

# Save the ground truth
save(groundTruth, file = output_file)
print(paste("Ground truth saved to:", output_file))