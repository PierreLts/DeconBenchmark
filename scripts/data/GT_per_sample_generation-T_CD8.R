#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop(paste("4 arguments must be supplied instead of", length(args)), call. = FALSE)
}

####### Parameter of script (ORDER IS IMPORTANT)
path_Rlibrary <- args[1] #IMPORTANT
seurat_file <- args[2]  # Path to Seurat RDS file
output_dir <- args[3]  # Output directory
prefix <- args[4]  # Prefix for output files

# Libraries
.libPaths(path_Rlibrary, FALSE) #IMPORTANT
library(Seurat)
library(reshape2)
library(dplyr)

# Load the Seurat object
print(paste("Loading Seurat object from:", seurat_file))
sc_data <- readRDS(seurat_file)

# Filter out "Undetermined" and "Multiplet" cells
valid_cells <- !sc_data@meta.data$Sample_Tag %in% c("Undetermined", "Multiplet")
filtered_sc <- subset(sc_data, cells = rownames(sc_data@meta.data)[valid_cells])

# Get the valid sample names
valid_sample_names <- unique(filtered_sc@meta.data$Sample_Name)
print(paste("Found", length(valid_sample_names), "valid samples"))

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
# Get all unique cell types across all samples
all_cell_types <- unique(unlist(lapply(sample_proportions, names)))
print(paste("Found", length(all_cell_types), "unique cell types"))

# Check if 'T_CD8_memory' and 'T_CD8_naive' exist
has_cd8_memory <- "T_CD8_memory" %in% all_cell_types
has_cd8_naive <- "T_CD8_naive" %in% all_cell_types
if (has_cd8_memory && has_cd8_naive) {
  print("Found both 'T_CD8_memory' and 'T_CD8_naive', will combine them")
} else if (has_cd8_memory || has_cd8_naive) {
  print("Warning: Found only one of 'T_CD8_memory' or 'T_CD8_naive', combining anyway")
} else {
  print("Warning: Neither 'T_CD8_memory' nor 'T_CD8_naive' found in cell types")
}

# Create a new cell type list with 'T_CD8' instead of 'T_CD8_memory' and 'T_CD8_naive'
new_cell_types <- all_cell_types
new_cell_types <- setdiff(new_cell_types, c("T_CD8_memory", "T_CD8_naive"))
new_cell_types <- c(new_cell_types, "T_CD8")  # Add the combined type
print(paste("New cell type list has", length(new_cell_types), "types"))

# Initialize an empty matrix
P <- matrix(0, 
            nrow = length(sample_proportions), 
            ncol = length(new_cell_types),
            dimnames = list(names(sample_proportions), new_cell_types))

# Fill in the matrix with the proportions
for (sample_name in names(sample_proportions)) {
  sample_props <- sample_proportions[[sample_name]]
  
  # For each cell type in the new list
  for (cell_type in new_cell_types) {
    if (cell_type == "T_CD8") {
      # Sum the values for 'T_CD8_memory' and 'T_CD8_naive'
      mem_val <- ifelse("T_CD8_memory" %in% names(sample_props), sample_props["T_CD8_memory"], 0)
      naive_val <- ifelse("T_CD8_naive" %in% names(sample_props), sample_props["T_CD8_naive"], 0)
      P[sample_name, cell_type] <- mem_val + naive_val
    } else if (cell_type %in% names(sample_props)) {
      # Copy existing values for other cell types
      P[sample_name, cell_type] <- sample_props[cell_type]
    }
    # Implicit else: cell type not present in this sample, keep as 0
  }
}

# Create a list to store the ground truth
groundTruth <- list(P = P)

# Print summary of combined matrix
print("Summary of ground truth proportions after combining T_CD8 subtypes:")
print(paste("Dimensions:", paste(dim(P), collapse=" x ")))
print("Cell types in the final matrix:")
print(colnames(P))

# Verify that all rows sum close to 1
row_sums <- rowSums(P)
print("Row sums (should be close to 1):")
print(range(row_sums))

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