#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
  stop("3 arguments must be supplied: R_LIBRARY_PATH OUTPUT_DIR PREFIX", call.=FALSE)
}

# Parameters
path_Rlibrary <- args[1]  # R library path
output_dir <- args[2]     # Output directory
prefix <- args[3]         # Prefix for file names

# Set library path
.libPaths(path_Rlibrary, FALSE)

# Print initial directory and permission info
print(paste("Working directory:", getwd()))
print(paste("Output directory:", output_dir))
print(paste("Directory exists:", dir.exists(output_dir)))
print(paste("Directory writable:", file.access(output_dir, 2) == 0))

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = TRUE, recursive = TRUE)

# Load libraries with error handling
tryCatch({
  library(DeconBenchmark)
  library(reshape2)
  print("CHECK: Libraries loaded successfully")
}, error = function(e) {
  print(paste("ERROR loading libraries:", e$message))
  quit(status = 1)
})

# Load BloodExample data with error handling
tryCatch({
  data(BloodExample)
  print("CHECK: BloodExample data loaded")
  
  # Verify data structure
  print(paste("BloodExample contains:", paste(names(BloodExample), collapse=", ")))
  print(paste("singleCellExpr dimensions:", paste(dim(BloodExample$singleCellExpr), collapse=" x ")))
  print(paste("singleCellLabels length:", length(BloodExample$singleCellLabels)))
}, error = function(e) {
  print(paste("ERROR loading BloodExample:", e$message))
  quit(status = 1)
})

print(paste("Using prefix:", prefix))

# Get the single cell data
singleCellExpr <- BloodExample$singleCellExpr
singleCellLabels <- BloodExample$singleCellLabels

# For BloodExample, we need to create simulated "samples" since it's a reference dataset
# We'll group cells by their cell type to create synthetic samples

# Get unique cell types
cell_types <- unique(singleCellLabels)
print(paste("Found", length(cell_types), "unique cell types"))
print(paste("Cell types:", paste(cell_types, collapse=", ")))

# Create a list to store the ground truth proportions for synthetic samples
sample_proportions <- list()

# For each cell type, create a "pure" sample (100% of that cell type)
for (ct in cell_types) {
  sample_name <- paste0("Pure_", ct)
  # Create a sample with 100% of this cell type
  props <- numeric(length(cell_types))
  names(props) <- cell_types
  props[ct] <- 1.0
  
  sample_proportions[[sample_name]] <- props
  print(paste("Created pure sample for", ct))
}

# Create a "Mixed" sample with equal proportion of all cell types
mixed_sample <- rep(1/length(cell_types), length(cell_types))
names(mixed_sample) <- cell_types
sample_proportions[["Mixed_Equal"]] <- mixed_sample
print("Created equally mixed sample")

# Create a "Weighted" mixed sample with proportions matching the original dataset
cell_counts <- table(singleCellLabels)
weighted_props <- as.numeric(cell_counts) / sum(cell_counts)
names(weighted_props) <- names(cell_counts)
sample_proportions[["Mixed_Weighted"]] <- weighted_props
print("Created weighted mixed sample")

# Verify the sample_proportions list
print(paste("Created", length(sample_proportions), "synthetic samples"))
print("Sample names:")
print(names(sample_proportions))

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

# Verify the matrix
print("Ground truth matrix dimensions:")
print(dim(P))
print("First few rows and columns of matrix P:")
print(P[1:min(3, nrow(P)), 1:min(3, ncol(P))])

# Create a list to store the ground truth
groundTruth <- list(P = P)

# Save as RDA with error handling
rda_filename <- file.path(output_dir, paste0(prefix, "_GT_proportions_per_sample.rda"))
tryCatch({
  save(groundTruth, file = rda_filename)
  print(paste("Successfully saved RDA to:", rda_filename))
  if (file.exists(rda_filename)) {
    print(paste("File exists with size:", file.info(rda_filename)$size, "bytes"))
  } else {
    print("WARNING: File doesn't exist after saving!")
  }
}, error = function(e) {
  print(paste("ERROR saving RDA file:", e$message))
})

# Save as CSV with error handling
csv_filename <- file.path(output_dir, paste0(prefix, "_GT_proportions_per_sample.csv"))
tryCatch({
  # Convert matrix to long format for CSV
  gt_df <- as.data.frame(groundTruth$P)
  gt_df$Sample <- rownames(gt_df)
  gt_long <- melt(gt_df, id.vars = "Sample", variable.name = "CellType", value.name = "Proportion")
  
  # Write CSV
  write.csv(gt_long, file = csv_filename, row.names = FALSE)
  print(paste("Successfully saved CSV to:", csv_filename))
  if (file.exists(csv_filename)) {
    print(paste("File exists with size:", file.info(csv_filename)$size, "bytes"))
  } else {
    print("WARNING: File doesn't exist after saving!")
  }
}, error = function(e) {
  print(paste("ERROR saving CSV file:", e$message))
})

# Final check of output directory
print("Files in output directory:")
print(list.files(output_dir, pattern = paste0(prefix, "_GT_proportions_per_sample")))

print(paste("Per-sample ground truth processing completed for prefix:", prefix))