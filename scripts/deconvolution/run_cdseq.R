#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop(paste("4 arguments must be supplied instead of", length(args)), call. = FALSE)
}

####### Parameter of script (ORDER IS IMPORTANT)
path_Rlibrary <- args[1] #IMPORTANT
input_dir <- args[2]     # Directory containing the generated data files
output_base_dir <- args[3]   # Base output directory for results
sample_filter <- args[4]     # Sample filter: A, B, or AB

# Libraries
.libPaths(path_Rlibrary, FALSE) #IMPORTANT
library(CDSeq)
library(DeconBenchmark)
library(Matrix)
print("CHECK: Libraries loaded")

# Extract dataset prefix from the input directory name
dataset_prefix <- basename(input_dir)

# Create dataset-specific output directory
output_dir <- file.path(output_base_dir, dataset_prefix)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
print(paste("Results will be saved to:", output_dir))
print(paste("Using sample filter:", sample_filter))

# Load individual files from the directory with filter suffix
bulk_path <- file.path(input_dir, paste0(dataset_prefix, "_bulk.rda"))
singleCellLabels_path <- file.path(input_dir, paste0(dataset_prefix, "_singleCellLabels_", sample_filter, ".rda"))
singleCellExpr_path <- file.path(input_dir, paste0(dataset_prefix, "_singleCellExpr_", sample_filter, ".rda"))

if (!file.exists(bulk_path)) {
  stop(paste("Bulk data file not found at:", bulk_path))
}
if (!file.exists(singleCellLabels_path)) {
  stop(paste("Single cell labels file not found at:", singleCellLabels_path))
}
if (!file.exists(singleCellExpr_path)) {
  stop(paste("SingleCellExpr file not found at:", singleCellExpr_path))
}

# Load the data
load(bulk_path)
load(singleCellLabels_path)
load(singleCellExpr_path)

print("CHECK: Data loaded successfully")

# Get common genes between bulk and single-cell data
common_genes <- intersect(rownames(bulk), rownames(singleCellExpr))
message(paste("Common genes between datasets:", length(common_genes)))
# Subset both matrices to common genes
bulk <- bulk[common_genes, , drop=FALSE]
singleCellExpr <- singleCellExpr[common_genes, , drop=FALSE]

# Get number of cell types from the single-cell data
n_cell_types <- length(unique(singleCellLabels))
print(paste("Number of cell types detected:", n_cell_types))

# Run CDSeq
print("Starting CDSeq deconvolution...")
cdseq.result <- CDSeq(
  bulk_data = bulk,
  cell_type_number = n_cell_types, #seq(max(2, n_cell_types-3), min(n_cell_types+3, 15), 3),
  beta = 0.5,                 
  alpha = 5,                           
  mcmc_iterations = 700,
  dilution_factor = 5,
  gene_subset_size = 300,
  block_number = 5,
  cpu_number = 16
)

print("CHECK: CDSeq deconvolution completed, now formatting...")

reference <- generateReference(singleCellExpr, singleCellLabels, c("markers", "sigGenes", "signature", "cellTypeExpr"), 1)
markers <- reference$markers
sigGenes <- reference$sigGenes
signature <- reference$signature
cellTypeExpr <- reference$cellTypeExpr

# For CDSeq cell type assignment, use the cellTypeExpr (avg expression by cell type)
print("Assigning cell types using DeconBenchmark reference profiles...")

# Calculate correlation matrix between estimated GEPs and reference profiles
corr_matrix <- cor(cdseq.result$estGEP, cellTypeExpr)

# Use cellTypeAssign with the correlation matrix
cell_assignment <- tryCatch({
  CDSeq::cellTypeAssign(corr_matrix, threshold = 0.7)
}, error = function(e) {
  message("Error with cellTypeAssign: ", e$message)
  
  # Create a simple manual mapping as fallback
  cell_types_idx <- apply(corr_matrix, 1, which.max)
  manual_mapping <- colnames(cellTypeExpr)[cell_types_idx]
  names(manual_mapping) <- colnames(cdseq.result$estGEP)
  return(list(mapping = manual_mapping))
})

# Structure the results
deconvolutionResult <- list()
deconvolutionResult$CDSeq <- list(
  P = t(cdseq.result$estProp),  # TRANSPOSED: Now samples in rows, cell types in columns
  S = cdseq.result$estGEP       # Genes in rows and cell types in columns
)

# Apply the cell type mapping
if (!is.null(cell_assignment) && !is.null(cell_assignment$mapping)) {
  # Use the mapping to rename columns
  new_colnames <- cell_assignment$mapping
  colnames(deconvolutionResult$CDSeq$P) <- new_colnames
  
  # Store the assignment info
  deconvolutionResult$CDSeq$cell_type_assignment <- cell_assignment
  print("Cell type assignment successful!")
} else {
  print("Cell type assignment failed or returned no mapping - using generic cell type names")
  colnames(deconvolutionResult$CDSeq$P) <- paste0("CellType", 1:ncol(deconvolutionResult$CDSeq$P))
}

# Normalize proportions if needed
prop_matrix <- deconvolutionResult$CDSeq$P

# Process each sample (row) individually
for (i in 1:nrow(prop_matrix)) {
  # Take absolute values to avoid negative proportions
  prop_matrix[i, ] <- abs(prop_matrix[i, ])
  
  # Normalize to sum to 1
  row_sum <- sum(prop_matrix[i, ])
  if (row_sum > 0) {  # Avoid division by zero
    prop_matrix[i, ] <- prop_matrix[i, ] / row_sum
  } else {
    # If all values are zero, keep as zeros
    prop_matrix[i, ] <- rep(0, ncol(prop_matrix))
    print(paste("Warning: Sample", rownames(prop_matrix)[i], "had all zero values - keeping as zeros"))
  }
}

# Update the proportions matrix in the result object
deconvolutionResult$CDSeq$P <- prop_matrix
print("Normalized cell type proportions (absolute values, sum to 1 per sample)")

# Save results as RDA with filter suffix
results_filename <- file.path(output_dir, paste0("results_CDSeq_", sample_filter, ".rda"))
save(deconvolutionResult, file=results_filename, compress=TRUE)

# Save proportions as CSV with filter suffix
proportions_df <- as.data.frame(deconvolutionResult$CDSeq$P)
proportions_df$Sample <- rownames(proportions_df)
csv_filename <- file.path(output_dir, paste0("results_CDSeq_", sample_filter, "_proportions.csv"))
write.csv(proportions_df, file=csv_filename, row.names=FALSE)

print(paste("Results saved to:", results_filename))
print(paste("Proportions CSV saved to:", csv_filename))