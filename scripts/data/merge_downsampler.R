#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 7) {
  stop("At least 7 arguments needed: R_LIBRARY_PATH OUTPUT_PREFIX CELLS_PER_TYPE FILTER_TYPE(A,B,AB) BATCH1_PREFIX BATCH2_PREFIX [BATCH3_PREFIX BATCH4_PREFIX]", call.=FALSE)
}

# Parameters
path_Rlibrary <- args[1]   # R library path
output_prefix <- args[2]   # Output prefix (e.g., "TB")
cells_per_type <- as.numeric(args[3])  # Number of cells per cell type
filter_type <- args[4]     # Filter type: "A", "B", or "AB"
batch_prefixes <- args[5:length(args)]  # Batch prefixes (e.g., "TB1", "TB2", etc.)

# Validate filter type
if (!filter_type %in% c("A", "B", "AB")) {
  stop("Filter type must be 'A', 'B', or 'AB'")
}

# Set library path
.libPaths(path_Rlibrary, FALSE)
library(Matrix)
library(dplyr)
library(reshape2)  # For melt function

# Calculate the output suffix based on batches and cells per type
num_batches <- length(batch_prefixes)
total_cells_per_type <- num_batches * cells_per_type
output_prefix_with_suffix <- paste0(output_prefix, "_D", total_cells_per_type)

# Base directories
base_dir <- "/work/gr-fe/lorthiois/DeconBenchmark/generated_data"
output_dir <- file.path(base_dir, output_prefix_with_suffix)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Build filter suffix based on filter type
filter_suffix <- paste0("_", filter_type)

cat("Will merge and downsample", length(batch_prefixes), "batches:", paste(batch_prefixes, collapse=", "), "\n")
cat("Cells per type per batch:", cells_per_type, "\n")
cat("Total cells per type:", total_cells_per_type, "\n")
cat("Using filter type:", filter_type, "\n")
cat("Output will be saved to:", output_dir, "\n")
cat("Output prefix:", output_prefix_with_suffix, "\n")

# === STEP 1: DOWN-SAMPLE AND MERGE LABELS AND SUBJECTS ===
cat("\n=== Downsampling and Merging Cell Labels and Subjects ===\n")

all_labels <- list()
all_subjects <- list()
all_full_labels <- list()  # For ground truth generation (without downsampling)
all_full_subjects <- list()

for (batch in batch_prefixes) {
  # Load labels
  labels_file <- file.path(base_dir, batch, paste0(batch, "_singleCellLabels", filter_suffix, ".rda"))
  if (file.exists(labels_file)) {
    load(labels_file)
    
    # Store full labels for ground truth generation - keep original cell names
    all_full_labels[[batch]] <- singleCellLabels
    
    # Down-sample the labels by cell type
    downsampled_labels <- c()
    for (cell_type in unique(singleCellLabels)) {
      # Get all cells of this type
      cells_of_type <- names(singleCellLabels[singleCellLabels == cell_type])
      # Sample cells or take all if fewer than requested
      if (length(cells_of_type) <= cells_per_type) {
        sampled_cells <- cells_of_type
        cat("  Warning: only", length(cells_of_type), "cells for type", cell_type, "in batch", batch, "(requested", cells_per_type, ")\n")
      } else {
        set.seed(42 + which(batch_prefixes == batch) + which(unique(singleCellLabels) == cell_type))  # Different seed for each batch and cell type
        sampled_cells <- sample(cells_of_type, cells_per_type)
      }
      # Add sampled cells to downsampled labels
      downsampled_labels[sampled_cells] <- singleCellLabels[sampled_cells]
    }
    
    all_labels[[batch]] <- downsampled_labels
    cat("Loaded and downsampled from", length(singleCellLabels), "to", length(downsampled_labels), "cell labels from", batch, "\n")
  } else {
    cat("WARNING: Labels file not found for", batch, " (", labels_file, ")\n")
  }
  
  # Load subjects
  subjects_file <- file.path(base_dir, batch, paste0(batch, "_singleCellSubjects", filter_suffix, ".rda"))
  if (file.exists(subjects_file)) {
    load(subjects_file)
    
    # Store full subjects for ground truth generation - keep original names
    all_full_subjects[[batch]] <- singleCellSubjects
    
    # Down-sample the subjects to match downsampled labels
    if (exists("downsampled_labels")) {
      downsampled_subjects <- singleCellSubjects[names(downsampled_labels)]
      all_subjects[[batch]] <- downsampled_subjects
      cat("Loaded and downsampled from", length(singleCellSubjects), "to", length(downsampled_subjects), "cell subjects from", batch, "\n")
    }
  } else {
    cat("WARNING: Subjects file not found for", batch, " (", subjects_file, ")\n")
  }
}

# Combine all downsampled labels
singleCellLabels <- unlist(all_labels, recursive = FALSE)
# Add batch information to subjects if present
singleCellSubjects <- unlist(all_subjects, recursive = FALSE)

# Also combine full labels and subjects for ground truth
full_singleCellLabels <- unlist(all_full_labels, recursive = FALSE)
full_singleCellSubjects <- unlist(all_full_subjects, recursive = FALSE)

# Report on merged data
cat("Merged downsampled labels count:", length(singleCellLabels), "\n")
cat("Unique cell types:", paste(unique(singleCellLabels), collapse=", "), "\n")
if (length(singleCellSubjects) > 0) {
  cat("Merged downsampled subjects count:", length(singleCellSubjects), "\n")
  cat("Unique subjects:", length(unique(singleCellSubjects)), "\n")
}

# Clean up cell names by removing any potential batch prefixes that might exist
clean_cell_names <- function(cell_names) {
  # This will remove any batchname prefixes like "TB1." or "TB2_" from cell names
  # It looks for patterns like "PREFIX." or "PREFIX_" at the start of cell names
  clean_names <- cell_names
  
  # Apply cleaning for each batch prefix
  for (batch in batch_prefixes) {
    # Check for various separator patterns (., _, -)
    clean_names <- gsub(paste0("^", batch, "\\."), "", clean_names)
    clean_names <- gsub(paste0("^", batch, "_"), "", clean_names)
    clean_names <- gsub(paste0("^", batch, "-"), "", clean_names)
  }
  
  return(clean_names)
}

# Get clean cell names for labels
clean_label_names <- clean_cell_names(names(singleCellLabels))
clean_subjects_names <- clean_cell_names(names(singleCellSubjects))

# Save merged labels with original cell names
save(singleCellLabels, file = file.path(output_dir, paste0(output_prefix_with_suffix, "_singleCellLabels", filter_suffix, ".rda")))

# Create data frame with clean cell names for CSV output
labels_df <- data.frame(
  CellBarcode = clean_label_names,
  CellType = singleCellLabels,
  stringsAsFactors = FALSE
)
write.csv(labels_df, file = file.path(output_dir, paste0(output_prefix_with_suffix, "_singleCellLabels", filter_suffix, ".csv")), 
          row.names = FALSE)

# Save merged subjects if available
if (length(singleCellSubjects) > 0) {
  save(singleCellSubjects, file = file.path(output_dir, paste0(output_prefix_with_suffix, "_singleCellSubjects", filter_suffix, ".rda")))
  
  # Create data frame with clean cell names for CSV output
  subjects_df <- data.frame(
    CellBarcode = clean_subjects_names,
    Subject = singleCellSubjects,
    stringsAsFactors = FALSE
  )
  write.csv(subjects_df, file = file.path(output_dir, paste0(output_prefix_with_suffix, "_singleCellSubjects", filter_suffix, ".csv")), 
            row.names = FALSE)
}

# === STEP 2: CREATE GROUND TRUTH FROM ALL LABELS (FULL, NOT DOWNSAMPLED) ===
cat("\n=== Creating Ground Truth from Full Labels (Not Downsampled) ===\n")

# Calculate overall proportions from full labels
cell_types <- unique(full_singleCellLabels)
cell_counts <- table(full_singleCellLabels)
proportions <- as.numeric(cell_counts) / sum(cell_counts)

# Create a matrix with cell types as columns (single row for overall average)
P <- matrix(proportions, nrow = 1, dimnames = list("true_proportions", names(cell_counts)))

# Create a list to store the ground truth
groundTruth <- list(P = P)

# Save ground truth
save(groundTruth, file = file.path(output_dir, paste0(output_prefix_with_suffix, "_GT_proportions.rda")))
# Save as CSV
gt_df <- as.data.frame(t(groundTruth$P))
gt_df$CellType <- rownames(gt_df)
colnames(gt_df) <- c("Proportion", "CellType")
gt_df <- gt_df[, c("CellType", "Proportion")]
write.csv(gt_df, file = file.path(output_dir, paste0(output_prefix_with_suffix, "_GT_proportions.csv")), 
          row.names = FALSE)

# === STEP 3: CALCULATE PER-SAMPLE GROUND TRUTH ===
cat("\n=== Creating Per-Sample Ground Truth (Based on Full Labels) ===\n")

# Create per-sample ground truth based on full subjects
if (length(full_singleCellSubjects) > 0) {
  # Group cells by subject
  subject_cells <- split(names(full_singleCellLabels), full_singleCellSubjects)
  
  # Initialize ground truth matrix
  sample_proportions <- matrix(0, 
                              nrow = length(subject_cells), 
                              ncol = length(cell_types),
                              dimnames = list(names(subject_cells), cell_types))
  
  # Calculate proportions for each sample
  for (subject in names(subject_cells)) {
    subject_cell_types <- full_singleCellLabels[subject_cells[[subject]]]
    counts <- table(factor(subject_cell_types, levels = cell_types))
    sample_proportions[subject, ] <- as.numeric(counts) / sum(counts)
  }
  
  # Save per-sample ground truth
  groundTruth <- list(P = sample_proportions)
  save(groundTruth, file = file.path(output_dir, paste0(output_prefix_with_suffix, "_GT_proportions_per_sample.rda")))
  
  # Save as CSV
  per_sample_df <- as.data.frame(groundTruth$P)
  per_sample_df$Sample <- rownames(per_sample_df)
  per_sample_long <- melt(per_sample_df, id.vars = "Sample", 
                                   variable.name = "CellType", value.name = "Proportion")
  write.csv(per_sample_long, file = file.path(output_dir, paste0(output_prefix_with_suffix, "_GT_proportions_per_sample.csv")), 
            row.names = FALSE)
  
  cat("Created per-sample ground truth for", nrow(sample_proportions), "samples\n")
}

# === STEP 4: MERGE AND DOWNSAMPLE EXPRESSION DATA ===
cat("\n=== Merging and Downsampling Expression Data ===\n")

# We need to create merged expression data for the downsampled cells
# First, we'll collect expression data from each batch for the downsampled cells

all_expr_data <- list()
common_genes <- NULL

for (batch in batch_prefixes) {
  # Get cell names for this batch
  batch_cells <- names(all_labels[[batch]])
  if (length(batch_cells) == 0) {
    cat("WARNING: No downsampled cells for batch", batch, "\n")
    next
  }
  
  # Load expression data
  expr_file <- file.path(base_dir, batch, paste0(batch, "_singleCellExpr", filter_suffix, ".rda"))
  if (file.exists(expr_file)) {
    load(expr_file)
    
    # Find cells that exist in both the expression data and our downsampled cells
    common_cells <- intersect(batch_cells, colnames(singleCellExpr))
    
    if (length(common_cells) > 0) {
      # Subset expression data to include only downsampled cells
      batch_expr <- singleCellExpr[, common_cells, drop=FALSE]
      
      # Store expression data with original column names
      all_expr_data[[batch]] <- batch_expr
      
      # Update common genes
      if (is.null(common_genes)) {
        common_genes <- rownames(batch_expr)
      } else {
        common_genes <- intersect(common_genes, rownames(batch_expr))
      }
      
      cat("Loaded expression data for", ncol(batch_expr), "cells from batch", batch, "\n")
    } else {
      cat("WARNING: No common cells found between downsampled labels and expression data for batch", batch, "\n")
    }
  } else {
    cat("WARNING: Expression file not found for", batch, " (", expr_file, ")\n")
  }
}

# Check if we have any expression data to merge
if (length(all_expr_data) > 0 && length(common_genes) > 0) {
  cat("Found", length(common_genes), "common genes across all batches\n")
  
  # Subset each batch's expression data to common genes
  for (batch in names(all_expr_data)) {
    all_expr_data[[batch]] <- all_expr_data[[batch]][common_genes, , drop=FALSE]
  }
  
  # Merge expression data from all batches
  # Using do.call cbind without modifying column names
  singleCellExpr <- do.call(cbind, all_expr_data)
  
  cat("Created merged expression matrix with dimensions:", nrow(singleCellExpr), "x", ncol(singleCellExpr), "\n")
  
  # Save the merged expression matrix
  save(singleCellExpr, file = file.path(output_dir, paste0(output_prefix_with_suffix, "_singleCellExpr", filter_suffix, ".rda")))
  
  # Save a small subset as CSV
  max_rows <- min(50, nrow(singleCellExpr))
  max_cols <- min(50, ncol(singleCellExpr))
  small_expr <- as.matrix(singleCellExpr[1:max_rows, 1:max_cols])
  write.csv(small_expr, file = file.path(output_dir, paste0(output_prefix_with_suffix, "_singleCellExpr", filter_suffix, "_50x50.csv")))
  
  cat("Saved merged expression data\n")
} else {
  cat("WARNING: Could not create merged expression data - insufficient data\n")
}

cat("\n=== Processing completed successfully ===\n")
cat("All merged and downsampled files saved to:", output_dir, "\n")
cat("Output prefix:", output_prefix_with_suffix, "\n")
cat("Filter type used:", filter_type, "\n")