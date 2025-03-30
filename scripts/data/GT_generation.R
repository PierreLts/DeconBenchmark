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
library(reshape2)  # For melt function

# Function to check if the cell types in a matrix need special handling
check_cell_types <- function(cell_types) {
  # Track which special cases are present
  info <- list(
    has_cd8_memory = "T_CD8_memory" %in% cell_types,
    has_cd8_naive = "T_CD8_naive" %in% cell_types,
    has_cd8 = "T_CD8" %in% cell_types,
    has_gamma_delta = "T_gamma_delta" %in% cell_types,
    has_other = "Other" %in% cell_types
  )
  
  # Print detected cell type patterns
  pattern_desc <- character(0)
  if(info$has_cd8_memory && info$has_cd8_naive && !info$has_cd8) {
    pattern_desc <- c(pattern_desc, "T_CD8 subtypes (memory/naive) without combined T_CD8")
  } else if(!info$has_cd8_memory && !info$has_cd8_naive && info$has_cd8) {
    pattern_desc <- c(pattern_desc, "Combined T_CD8 without subtypes")
  } else if(info$has_cd8_memory && info$has_cd8_naive && info$has_cd8) {
    pattern_desc <- c(pattern_desc, "Both T_CD8 subtypes AND combined T_CD8 (unusual)")
  }
  
  if(info$has_gamma_delta && !info$has_other) {
    pattern_desc <- c(pattern_desc, "T_gamma_delta without Other")
  } else if(!info$has_gamma_delta && info$has_other) {
    pattern_desc <- c(pattern_desc, "Other without T_gamma_delta")
  } else if(info$has_gamma_delta && info$has_other) {
    pattern_desc <- c(pattern_desc, "Both T_gamma_delta AND Other")
  }
  
  if(length(pattern_desc) > 0) {
    print(paste("Cell type pattern detected:", paste(pattern_desc, collapse=", ")))
  }
  
  return(info)
}

# Load per-sample ground truth file - require this file
per_sample_gt_path <- file.path(input_dir, paste0(prefix, "_GT_proportions_per_sample.rda"))
if (!file.exists(per_sample_gt_path)) {
  stop(paste("Error: Per-sample ground truth file not found:", per_sample_gt_path))
}

print(paste("Loading per-sample ground truth from:", per_sample_gt_path))
load(per_sample_gt_path)

# Verify groundTruth object exists with the required structure
if (!exists("groundTruth") || is.null(groundTruth$P)) {
  stop("Error: Expected 'groundTruth$P' object not found in per-sample file")
}

# This is the per-sample matrix
per_sample_proportions <- groundTruth$P

# Print information about the per-sample data
print("Per-sample ground truth summary:")
print(paste("Samples:", nrow(per_sample_proportions)))
print(paste("Cell types:", paste(colnames(per_sample_proportions), collapse=", ")))

# Check cell type patterns
cell_type_info <- check_cell_types(colnames(per_sample_proportions))

# Calculate overall proportions by averaging across all samples
# This gives equal weight to each sample regardless of cell count
overall_proportions <- colMeans(per_sample_proportions)

# Create a new matrix for overall proportions
P <- matrix(overall_proportions, nrow = 1, 
            dimnames = list("true_proportions", names(overall_proportions)))

# Create a list to store the ground truth
groundTruth <- list(P = P)

# Print information for verification
print("Ground truth proportions derived from per-sample data:")
print(groundTruth$P)
print(paste("Total cell types:", length(colnames(P))))
print(paste("Sum of proportions:", sum(overall_proportions)))

# Save as RDA
rda_filename <- file.path(output_dir, paste0(prefix, "_GT_proportions.rda"))
save(groundTruth, file = rda_filename)

# Save as CSV
csv_filename <- file.path(output_dir, paste0(prefix, "_GT_proportions.csv"))
gt_df <- as.data.frame(t(groundTruth$P))
gt_df$CellType <- rownames(gt_df)
colnames(gt_df) <- c("Proportion", "CellType")
gt_df <- gt_df[, c("CellType", "Proportion")]
write.csv(gt_df, file = csv_filename, row.names = FALSE)

print(paste("Ground truth saved to:", rda_filename, "and", csv_filename))