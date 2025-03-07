#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop(paste("4 arguments must be supplied instead of", length(args)), call. = FALSE)
}

####### Parameter of script (ORDER IS IMPORTANT)
path_Rlibrary <- args[1] #IMPORTANT
input_dir <- args[2]     # Directory containing the input data
output_dir <- args[3]    # Output directory
prefix <- args[4]        # Prefix for output files

# Libraries
.libPaths(path_Rlibrary, FALSE) #IMPORTANT

# Load single cell labels file
sc_labels_path <- file.path(input_dir, paste0(prefix, "_singleCellLabels.rda"))
if (!file.exists(sc_labels_path)) {
  stop(paste("Single cell labels file not found:", sc_labels_path))
}

# Load the single cell labels
load(sc_labels_path)

# Calculate the ground truth proportions from the labels
if (!exists("singleCellLabels")) {
  stop("singleCellLabels object not found in the loaded file")
}

# Calculate overall proportions
cell_types <- unique(singleCellLabels)
cell_counts <- table(singleCellLabels)
proportions <- as.numeric(cell_counts) / sum(cell_counts)

# Create a matrix with cell types as columns (single row for overall average)
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