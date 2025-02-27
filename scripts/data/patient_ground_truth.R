#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
  stop("3 arguments must be supplied: <R_library_path> <input_scRNAseq_data> <output_dir>", call.=FALSE)
}

####### Parameters
path_Rlibrary <- args[1]
input_data <- args[2]
output_dir <- args[3]

# Libraries
.libPaths(path_Rlibrary, FALSE)
library(Seurat)
library(dplyr)
library(Matrix)

# Load Seurat object
print("Loading Seurat object...")
seurat_obj <- readRDS(input_data)

# Filter out undetermined and multiplet cells
print("Filtering cells...")
cell_filter <- !seurat_obj@meta.data$Sample_Tag %in% c("Undetermined", "Multiplet")
filtered_obj <- subset(seurat_obj, cells = rownames(seurat_obj@meta.data)[cell_filter])

# Extract metadata we need
meta_data <- filtered_obj@meta.data
cell_types <- meta_data$Cell_Type_Experimental
sample_ids <- meta_data$Sample_Tag  # This should contain patient identifiers

# Get unique patient IDs and cell types
unique_patients <- unique(sample_ids)
unique_cell_types <- unique(cell_types)

print(paste("Found", length(unique_patients), "unique patients"))
print(paste("Found", length(unique_cell_types), "unique cell types"))

# Create empty matrix to store proportions
# Rows = patients, Columns = cell types
ground_truth_matrix <- matrix(0, 
                              nrow = length(unique_patients), 
                              ncol = length(unique_cell_types))
rownames(ground_truth_matrix) <- unique_patients
colnames(ground_truth_matrix) <- unique_cell_types

# Calculate proportions for each patient
print("Calculating cell type proportions per patient...")
for (patient in unique_patients) {
  # Get cells for this patient
  patient_cells <- rownames(meta_data)[sample_ids == patient]
  total_cells <- length(patient_cells)
  
  # Skip if no cells found (shouldn't happen after filtering)
  if (total_cells == 0) {
    warning(paste("No cells found for patient", patient))
    next
  }
  
  # Get cell types for this patient
  patient_cell_types <- cell_types[sample_ids == patient]
  
  # Calculate proportions
  for (cell_type in unique_cell_types) {
    count <- sum(patient_cell_types == cell_type)
    proportion <- count / total_cells
    ground_truth_matrix[patient, cell_type] <- proportion
  }
  
  # Verify proportions sum to 1
  row_sum <- sum(ground_truth_matrix[patient, ])
  if (abs(row_sum - 1) > 1e-10) {
    warning(paste("Proportions for patient", patient, "sum to", row_sum, "not 1"))
  }
}

# Save results
print("Saving patient-specific ground truth matrix...")
groundTruth <- list(
  P = ground_truth_matrix,
  patients = unique_patients,
  cell_types = unique_cell_types
)

# Create feature mapping
# Extract patient IDs from Sample_Tag (assuming format like "01A" where "01" is patient ID)
patient_mapping <- data.frame(
  sample = rownames(ground_truth_matrix),
  patient_id = gsub("[AB]$", "", rownames(ground_truth_matrix)),
  time = ifelse(grepl("A$", rownames(ground_truth_matrix)), "before", "after"),
  stringsAsFactors = FALSE
)

# Add mapping to groundTruth
groundTruth$patient_mapping <- patient_mapping

# Save to output directory
output_file <- file.path(output_dir, "patient_ground_truth_proportions.rda")
save(groundTruth, file = output_file)
print(paste("Ground truth proportions saved to:", output_file))

# Also save a CSV version for easier inspection
write.csv(ground_truth_matrix, file = file.path(output_dir, "patient_ground_truth_proportions.csv"))
print("Done!")