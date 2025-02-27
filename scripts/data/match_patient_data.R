#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop("4 arguments must be supplied: <R_library_path> <features_csv> <ground_truth_file> <output_dir>", call.=FALSE)
}

####### Parameters
path_Rlibrary <- args[1]
features_file <- args[2]
ground_truth_file <- args[3]
output_dir <- args[4]

# Libraries
.libPaths(path_Rlibrary, FALSE)
library(dplyr)

# Load data
features <- read.csv(features_file, stringsAsFactors = FALSE)
load(ground_truth_file)  # Should load 'groundTruth' object

# Create mapping between bulk samples and patient identifiers
print("Creating bulk-to-scRNA sample mapping...")

# Extract bulk sample IDs and match with ground truth patient IDs
patient_mapping <- groundTruth$patient_mapping

# Match with features data
matched_data <- merge(
  features,
  patient_mapping,
  by = c("sample", "time"),
  all.x = FALSE,  # Keep only matched samples
  all.y = FALSE   # Keep only matched patients
)

# Save the mapping
matched_file <- file.path(output_dir, "patient_sample_mapping.csv")
write.csv(matched_data, file = matched_file, row.names = FALSE)
print(paste("Patient-sample mapping saved to:", matched_file))

# Create patient-specific ground truth matrices for later use in benchmarking
# This will create one RDA file per patient that contains their ground truth proportions

for (pid in unique(matched_data$patient_id)) {
  # Get samples for this patient
  patient_samples <- matched_data$sample[matched_data$patient_id == pid]
  
  # Get ground truth rows for this patient
  patient_gt_indices <- which(gsub("[AB]$", "", rownames(groundTruth$P)) == pid)
  
  if (length(patient_gt_indices) == 0) {
    warning(paste("No ground truth found for patient", pid))
    next
  }
  
  # Extract ground truth for this patient
  patient_ground_truth <- groundTruth$P[patient_gt_indices, , drop = FALSE]
  
  # Create patient-specific ground truth object
  patientGroundTruth <- list(
    P = patient_ground_truth,
    patient_id = pid,
    samples = patient_samples,
    cell_types = colnames(patient_ground_truth)
  )
  
  # Save to output directory
  output_file <- file.path(output_dir, paste0("ground_truth_patient_", pid, ".rda"))
  save(patientGroundTruth, file = output_file)
  print(paste("Ground truth for patient", pid, "saved to:", output_file))
}

# Also create a combined ground truth file with matched samples only
matched_gt <- groundTruth$P[rownames(groundTruth$P) %in% matched_data$sample, , drop = FALSE]
matchedGroundTruth <- list(
  P = matched_gt,
  samples = rownames(matched_gt),
  cell_types = colnames(matched_gt),
  patient_mapping = matched_data
)

# Save matched ground truth
matched_output_file <- file.path(output_dir, "matched_ground_truth_proportions.rda")
save(matchedGroundTruth, file = matched_output_file)
print(paste("Matched ground truth saved to:", matched_output_file))

print("Done!")