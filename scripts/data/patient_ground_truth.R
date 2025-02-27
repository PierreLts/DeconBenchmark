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
sample_ids <- meta_data$Sample_Tag  # Contains patient identifiers (e.g., "01A", "02B")

# Get unique sample IDs and cell types
unique_samples <- unique(sample_ids)
unique_cell_types <- unique(cell_types)

# Extract patient IDs and timepoints from sample IDs
extract_patient_info <- function(sample_id) {
  # Assuming format like "01A" where digits are patient ID and letter is timepoint
  patient_id <- gsub("[AB]$", "", sample_id)
  timepoint <- ifelse(grepl("A$", sample_id), "before", "after")
  return(list(patient_id = patient_id, timepoint = timepoint))
}

# Process sample IDs to extract patient information
patient_info <- lapply(unique_samples, extract_patient_info)
names(patient_info) <- unique_samples

# Create a more structured mapping
patient_mapping <- data.frame(
  sample = unique_samples,
  patient_id = sapply(patient_info, function(x) x$patient_id),
  timepoint = sapply(patient_info, function(x) x$timepoint),
  stringsAsFactors = FALSE
)

print(paste("Found", length(unique(patient_mapping$patient_id)), "unique patients"))
print(paste("Found", length(unique_cell_types), "unique cell types"))

# Calculate cell type proportions for each sample
print("Calculating cell type proportions per sample...")
ground_truth_matrix <- matrix(0, 
                              nrow = length(unique_samples), 
                              ncol = length(unique_cell_types))
rownames(ground_truth_matrix) <- unique_samples
colnames(ground_truth_matrix) <- unique_cell_types

for (sample_id in unique_samples) {
  # Get cells for this sample
  sample_cells <- rownames(meta_data)[sample_ids == sample_id]
  total_cells <- length(sample_cells)
  
  # Skip if no cells found
  if (total_cells == 0) {
    warning(paste("No cells found for sample", sample_id))
    next
  }
  
  # Get cell types for this sample
  sample_cell_types <- cell_types[sample_ids == sample_id]
  
  # Calculate proportions
  for (cell_type in unique_cell_types) {
    count <- sum(sample_cell_types == cell_type)
    proportion <- count / total_cells
    ground_truth_matrix[sample_id, cell_type] <- proportion
  }
  
  # Verify proportions sum to 1
  row_sum <- sum(ground_truth_matrix[sample_id, ])
  if (abs(row_sum - 1) > 1e-10) {
    warning(paste("Proportions for sample", sample_id, "sum to", row_sum, "not 1"))
  }
}

# Save the ground truth proportions and patient mapping
groundTruth <- list(
  P = ground_truth_matrix,
  patient_mapping = patient_mapping,
  cell_types = unique_cell_types
)

# Save to output directory
output_file <- file.path(output_dir, "patient_ground_truth_proportions.rda")
save(groundTruth, file = output_file)
print(paste("Ground truth proportions saved to:", output_file))

# Also generate per-patient ground truth files
unique_patients <- unique(patient_mapping$patient_id)
for (patient_id in unique_patients) {
  # Get samples for this patient
  patient_samples <- patient_mapping$sample[patient_mapping$patient_id == patient_id]
  
  # Extract ground truth for these samples
  patient_gt <- ground_truth_matrix[patient_samples, , drop = FALSE]
  
  # Create patient-specific ground truth object
  patientGroundTruth <- list(
    P = patient_gt,
    patient_id = patient_id,
    sample_ids = patient_samples,
    cell_types = unique_cell_types,
    timepoints = patient_mapping$timepoint[patient_mapping$patient_id == patient_id]
  )
  
  # Save to output directory
  patient_file <- file.path(output_dir, paste0("ground_truth_patient_", patient_id, ".rda"))
  save(patientGroundTruth, file = patient_file)
  print(paste("Ground truth for patient", patient_id, "saved to:", patient_file))
}

# Also save a CSV version for easier inspection
write.csv(ground_truth_matrix, file = file.path(output_dir, "patient_ground_truth_proportions.csv"))
print("Done!")