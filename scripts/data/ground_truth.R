#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
  stop(paste("3 arguments must be supplied instead of", length(args)), call. = FALSE)
}

####### Parameter of script (ORDER IS IMPORTANT)
path_Rlibrary <- args[1] #IMPORTANT
input_dir <- args[2]  # Directory containing the generated data
output_file <- args[3]  # Full path to output ground truth file

# Libraries
.libPaths(path_Rlibrary, FALSE) #IMPORTANT

# Find and load the Batch1.rda file which should contain all necessary data
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
  # Load the batch file which should contain all data
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

# Save the ground truth
save(groundTruth, file = output_file)
print(paste("Ground truth saved to:", output_file))