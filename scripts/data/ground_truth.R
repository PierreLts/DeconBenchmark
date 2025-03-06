#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop(paste("4 arguments must be supplied instead of", length(args)), call. = FALSE)
}

####### Parameter of script (ORDER IS IMPORTANT)
path_Rlibrary <- args[1] #IMPORTANT
input_dir <- args[2]  # Directory containing the generated data
output_dir <- args[3]  # Output directory
prefix <- args[4]  # Prefix for output files

# Libraries
.libPaths(path_Rlibrary, FALSE) #IMPORTANT

# Load Batch1.rda
batch_file <- file.path(input_dir, "Batch1.rda")
if (!file.exists(batch_file)) {
  # Try to find alternative data files
  data_files <- list.files(input_dir, pattern = "\\.rda$", full.names = TRUE)
  sc_labels_file <- grep("singleCellLabels", data_files, value = TRUE)
  
  if (length(sc_labels_file) > 0) {
    # Load single cell labels directly
    load(sc_labels_file[1])
  } else {
    stop("Cannot find Batch1.rda or singleCellLabels.rda in the input directory")
  }
} else {
  load(batch_file)
  
  # Extract single cell labels from the Batch1 object
  if (exists("Batch1") && is.list(Batch1) && !is.null(Batch1$singleCellLabels)) {
    singleCellLabels <- Batch1$singleCellLabels
  } else {
    stop("Batch1 object does not contain singleCellLabels")
  }
}

# Check if singleCellLabels exists
if (!exists("singleCellLabels")) {
  stop("singleCellLabels not found in the input data")
}

# Calculate the ground truth proportions by counting cell types
cell_types <- unique(singleCellLabels)
cell_counts <- table(singleCellLabels)
proportions <- as.numeric(cell_counts) / sum(cell_counts)

# Create a matrix with cell types as columns
P <- matrix(proportions, nrow = 1, dimnames = list("true_proportions", names(cell_counts)))

# Create a list to store the ground truth
groundTruth <- list(P = P)

# Print information for verification
print("Ground truth proportions:")
print(groundTruth$P)
print(paste("Total cell types:", length(cell_types)))
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