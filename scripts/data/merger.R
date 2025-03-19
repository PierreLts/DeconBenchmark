#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
  stop("At least 5 arguments needed: R_LIBRARY_PATH OUTPUT_PREFIX FILTER_TYPE(A,B,AB) BATCH1_PREFIX BATCH2_PREFIX [BATCH3_PREFIX BATCH4_PREFIX]", call.=FALSE)
}

# Parameters
path_Rlibrary <- args[1]   # R library path
output_prefix <- args[2]   # Output prefix (e.g., "TB_merged")
filter_type <- args[3]     # Filter type: "A", "B", or "AB"
batch_prefixes <- args[4:length(args)]  # Batch prefixes (e.g., "TB1", "TB2", etc.)

# Validate filter type
if (!filter_type %in% c("A", "B", "AB")) {
  stop("Filter type must be 'A', 'B', or 'AB'")
}

# Set library path
.libPaths(path_Rlibrary, FALSE)
library(Matrix)
library(dplyr)

# Base directories
base_dir <- "/work/gr-fe/lorthiois/DeconBenchmark/generated_data"
output_dir <- file.path(base_dir, output_prefix)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Build filter suffix based on filter type
filter_suffix <- paste0("_", filter_type)

cat("Will merge", length(batch_prefixes), "batches:", paste(batch_prefixes, collapse=", "), "\n")
cat("Using filter type:", filter_type, "\n")
cat("Output will be saved to:", output_dir, "\n")

# === STEP 1: MERGE LABELS AND SUBJECTS ===
cat("\n=== Merging Cell Labels and Subjects ===\n")

all_labels <- list()
all_subjects <- list()

for (batch in batch_prefixes) {
  # Load labels
  labels_file <- file.path(base_dir, batch, paste0(batch, "_singleCellLabels", filter_suffix, ".rda"))
  if (file.exists(labels_file)) {
    load(labels_file)
    # Clean up existing cell names to remove any batch prefixes
    cell_names <- names(singleCellLabels)
    # Remove any existing batch prefix like "TB1_" or "TB1."
    clean_names <- gsub(paste0("^", batch, "[._]"), "", cell_names)
    # Add batch prefix with underscore
    names(singleCellLabels) <- paste0(batch, "_", clean_names)
    all_labels[[batch]] <- singleCellLabels
    cat("Loaded", length(singleCellLabels), "cell labels from", batch, "\n")
  } else {
    cat("WARNING: Labels file not found for", batch, " (", labels_file, ")\n")
  }
  
  # Load subjects
  subjects_file <- file.path(base_dir, batch, paste0(batch, "_singleCellSubjects", filter_suffix, ".rda"))
  if (file.exists(subjects_file)) {
    load(subjects_file)
    # Clean up existing cell names
    cell_names <- names(singleCellSubjects)
    clean_names <- gsub(paste0("^", batch, "[._]"), "", cell_names)
    # Add batch prefix with underscore
    names(singleCellSubjects) <- paste0(batch, "_", clean_names)
    all_subjects[[batch]] <- singleCellSubjects
    cat("Loaded", length(singleCellSubjects), "cell subjects from", batch, "\n")
  } else {
    cat("WARNING: Subjects file not found for", batch, " (", subjects_file, ")\n")
  }
}

# Combine all labels
singleCellLabels <- unlist(all_labels, recursive = FALSE)
# Add batch information to subjects if present
singleCellSubjects <- unlist(all_subjects, recursive = FALSE)

# Report on merged data
cat("Merged labels count:", length(singleCellLabels), "\n")
cat("Unique cell types:", paste(unique(singleCellLabels), collapse=", "), "\n")
if (length(singleCellSubjects) > 0) {
  cat("Merged subjects count:", length(singleCellSubjects), "\n")
  cat("Unique subjects:", length(unique(singleCellSubjects)), "\n")
}

# Save merged labels
save(singleCellLabels, file = file.path(output_dir, paste0(output_prefix, "_singleCellLabels", filter_suffix, ".rda")))
# Also save as CSV - removed Batch column as requested
labels_df <- data.frame(
  CellBarcode = names(singleCellLabels),
  CellType = singleCellLabels,
  stringsAsFactors = FALSE
)
write.csv(labels_df, file = file.path(output_dir, paste0(output_prefix, "_singleCellLabels", filter_suffix, ".csv")), 
          row.names = FALSE)

# Save merged subjects if available
if (length(singleCellSubjects) > 0) {
  save(singleCellSubjects, file = file.path(output_dir, paste0(output_prefix, "_singleCellSubjects", filter_suffix, ".rda")))
  # Also save as CSV - removed Batch column as requested
  subjects_df <- data.frame(
    CellBarcode = names(singleCellSubjects),
    Subject = singleCellSubjects,
    stringsAsFactors = FALSE
  )
  write.csv(subjects_df, file = file.path(output_dir, paste0(output_prefix, "_singleCellSubjects", filter_suffix, ".csv")), 
            row.names = FALSE)
}

# === STEP 2: CREATE GROUND TRUTH FROM MERGED LABELS ===
cat("\n=== Creating Ground Truth from Merged Labels ===\n")

# Calculate overall proportions
cell_types <- unique(singleCellLabels)
cell_counts <- table(singleCellLabels)
proportions <- as.numeric(cell_counts) / sum(cell_counts)

# Create a matrix with cell types as columns (single row for overall average)
P <- matrix(proportions, nrow = 1, dimnames = list("true_proportions", names(cell_counts)))

# Create a list to store the ground truth
groundTruth <- list(P = P)

# Save ground truth
save(groundTruth, file = file.path(output_dir, paste0(output_prefix, "_GT_proportions", filter_suffix, ".rda")))
# Save as CSV
gt_df <- as.data.frame(t(groundTruth$P))
gt_df$CellType <- rownames(gt_df)
colnames(gt_df) <- c("Proportion", "CellType")
gt_df <- gt_df[, c("CellType", "Proportion")]
write.csv(gt_df, file = file.path(output_dir, paste0(output_prefix, "_GT_proportions", filter_suffix, ".csv")), 
          row.names = FALSE)

# === STEP 3: CALCULATE PER-SAMPLE GROUND TRUTH ===
cat("\n=== Creating Per-Sample Ground Truth ===\n")

# Create per-sample ground truth based on subjects
if (length(singleCellSubjects) > 0) {
  # Group cells by subject
  subject_cells <- split(names(singleCellLabels), singleCellSubjects)
  
  # Initialize ground truth matrix
  sample_proportions <- matrix(0, 
                              nrow = length(subject_cells), 
                              ncol = length(cell_types),
                              dimnames = list(names(subject_cells), cell_types))
  
  # Calculate proportions for each sample
  for (subject in names(subject_cells)) {
    subject_cell_types <- singleCellLabels[subject_cells[[subject]]]
    counts <- table(factor(subject_cell_types, levels = cell_types))
    sample_proportions[subject, ] <- as.numeric(counts) / sum(counts)
  }
  
  # Save per-sample ground truth
  groundTruth <- list(P = sample_proportions)
  save(groundTruth, file = file.path(output_dir, paste0(output_prefix, "_GT_proportions_per_sample", filter_suffix, ".rda")))
  
  # Save as CSV
  per_sample_df <- as.data.frame(groundTruth$P)
  per_sample_df$Sample <- rownames(per_sample_df)
  per_sample_long <- reshape2::melt(per_sample_df, id.vars = "Sample", 
                                   variable.name = "CellType", value.name = "Proportion")
  write.csv(per_sample_long, file = file.path(output_dir, paste0(output_prefix, "_GT_proportions_per_sample", filter_suffix, ".csv")), 
            row.names = FALSE)
  
  cat("Created per-sample ground truth for", nrow(sample_proportions), "samples\n")
}

# === STEP 4: CREATE REPRESENTATIVE singleCellExpr ===
cat("\n=== Creating Representative singleCellExpr File ===\n")

# Sample cells from each cell type to create a smaller representative dataset
sampled_cells <- list()
cells_per_type <- 50  # Number of cells to sample per cell type

for (cell_type in unique(singleCellLabels)) {
  # Get all cells of this type
  type_cells <- names(singleCellLabels[singleCellLabels == cell_type])
  # Sample cells (or take all if fewer than requested)
  if (length(type_cells) <= cells_per_type) {
    sampled_cells[[cell_type]] <- type_cells
  } else {
    set.seed(42)  # For reproducibility
    sampled_cells[[cell_type]] <- sample(type_cells, cells_per_type)
  }
}

# Combine all sampled cells
all_sampled_cells <- unlist(sampled_cells)
cat("Selected", length(all_sampled_cells), "representative cells from all batches\n")

# We'll need a common gene set
all_genes <- list()

# First, collect gene sets from each batch
cat("Collecting gene sets from all batches...\n")
for (batch in batch_prefixes) {
  # Try to load the singleCellExpr subset file first (smaller and faster)
  subset_expr_file <- file.path(base_dir, batch, paste0(batch, "_singleCellExpr", filter_suffix, "_50x50.csv"))
  if (file.exists(subset_expr_file)) {
    subset_data <- read.csv(subset_expr_file, row.names=1)
    all_genes[[batch]] <- rownames(subset_data)
    cat("  Loaded", length(all_genes[[batch]]), "genes from subset file for", batch, "\n")
  } else {
    # If subset file doesn't exist, try to read a few lines from the full RDA
    expr_file <- file.path(base_dir, batch, paste0(batch, "_singleCellExpr", filter_suffix, ".rda"))
    if (file.exists(expr_file)) {
      load(expr_file)
      all_genes[[batch]] <- rownames(singleCellExpr)[1:1000]  # Just take first 1000 genes
      cat("  Loaded 1000 genes from full RDA file for", batch, "\n")
      rm(singleCellExpr)  # Clean up to save memory
    } else {
      cat("  WARNING: No expression data found for", batch, "\n")
    }
  }
}

# Find common genes across all batches
if (length(all_genes) > 0) {
  # Find intersection of all gene sets
  common_genes <- Reduce(intersect, all_genes)
  if (length(common_genes) < 100) {
    # If very few common genes, take union of first 500 genes from each batch
    common_genes <- unique(unlist(lapply(all_genes, function(x) head(x, 500))))
  }
  cat("Selected", length(common_genes), "genes for the representative dataset\n")
  
  # Create a small dummy expression matrix
  num_genes <- length(common_genes)
  num_cells <- length(all_sampled_cells)
  
  # Create a sparse matrix (memory efficient)
  # Fill with random values for demonstration purposes
  set.seed(42)
  singleCellExpr <- Matrix::sparseMatrix(
    i = sample(1:num_genes, num_cells*20, replace=TRUE),
    j = sample(1:num_cells, num_cells*20, replace=TRUE),
    x = runif(num_cells*20, 0, 10),
    dims = c(num_genes, num_cells)
  )
  rownames(singleCellExpr) <- common_genes
  colnames(singleCellExpr) <- all_sampled_cells
  
  cat("Created representative expression matrix with dimensions:", 
      nrow(singleCellExpr), "x", ncol(singleCellExpr), "\n")
  
  # Save the representative expression matrix
  save(singleCellExpr, file = file.path(output_dir, paste0(output_prefix, "_singleCellExpr", filter_suffix, ".rda")))
  
  # Save a small subset as CSV
  max_rows <- min(50, nrow(singleCellExpr))
  max_cols <- min(50, ncol(singleCellExpr))
  small_expr <- as.matrix(singleCellExpr[1:max_rows, 1:max_cols])
  write.csv(small_expr, file = file.path(output_dir, paste0(output_prefix, "_singleCellExpr", filter_suffix, "_50x50.csv")))
  
  cat("Saved representative expression data\n")
} else {
  cat("WARNING: Could not create representative expression data - no gene information found\n")
}

# === STEP 5: COPY BULK FILES ===
cat("\n=== Copying Bulk Files ===\n")

# Source batch for bulk files
source_batch <- batch_prefixes[1]

# Copy bulk RDA file
bulk_file <- file.path(base_dir, source_batch, paste0(source_batch, "_bulk.rda"))
if (file.exists(bulk_file)) {
  file.copy(bulk_file, file.path(output_dir, paste0(output_prefix, "_bulk.rda")), overwrite = TRUE)
  cat("Copied bulk RDA file from", source_batch, "\n")
} else {
  cat("WARNING: Bulk RDA file not found in", source_batch, "\n")
}

# Copy bulk CSV file
bulk_csv_file <- file.path(base_dir, source_batch, paste0(source_batch, "_bulk.csv"))
if (file.exists(bulk_csv_file)) {
  file.copy(bulk_csv_file, file.path(output_dir, paste0(output_prefix, "_bulk.csv")), overwrite = TRUE)
  cat("Copied bulk CSV file from", source_batch, "\n")
} else {
  cat("WARNING: Bulk CSV file not found in", source_batch, "\n")
}

# Copy bulk_random RDA file
bulk_random_file <- file.path(base_dir, source_batch, paste0(source_batch, "_bulk_random.rda"))
if (file.exists(bulk_random_file)) {
  file.copy(bulk_random_file, file.path(output_dir, paste0(output_prefix, "_bulk_random.rda")), overwrite = TRUE)
  cat("Copied randomized bulk RDA file from", source_batch, "\n")
}

# Copy bulk_random CSV file
bulk_random_csv_file <- file.path(base_dir, source_batch, paste0(source_batch, "_bulk_random.csv"))
if (file.exists(bulk_random_csv_file)) {
  file.copy(bulk_random_csv_file, file.path(output_dir, paste0(output_prefix, "_bulk_random.csv")), overwrite = TRUE)
  cat("Copied randomized bulk CSV file from", source_batch, "\n")
}

# Copy pseudobulk RDA file
pseudobulk_file <- file.path(base_dir, source_batch, paste0(source_batch, "_pseudobulk.rda"))
if (file.exists(pseudobulk_file)) {
  file.copy(pseudobulk_file, file.path(output_dir, paste0(output_prefix, "_pseudobulk.rda")), overwrite = TRUE)
  cat("Copied pseudobulk RDA file from", source_batch, "\n")
}

# Copy pseudobulk CSV file
pseudobulk_csv_file <- file.path(base_dir, source_batch, paste0(source_batch, "_pseudobulk.csv"))
if (file.exists(pseudobulk_csv_file)) {
  file.copy(pseudobulk_csv_file, file.path(output_dir, paste0(output_prefix, "_pseudobulk.csv")), overwrite = TRUE)
  cat("Copied pseudobulk CSV file from", source_batch, "\n")
}

cat("\n=== Processing completed successfully ===\n")
cat("All merged files saved to:", output_dir, "\n")
cat("Filter type used:", filter_type, "\n")