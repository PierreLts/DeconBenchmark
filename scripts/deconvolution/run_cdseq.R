#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
  stop(paste("5 arguments must be supplied instead of", length(args)), call. = FALSE)
}

####### Parameter of script (ORDER IS IMPORTANT)
path_Rlibrary <- args[1] #IMPORTANT
input_dir <- args[2]     # Directory containing the generated data files
output_base_dir <- args[3]   # Base output directory for results
sample_filter <- args[4]     # Sample filter: A, B, or AB
bulk_type <- args[5]         # New parameter: bulk file type

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
print(paste("Using bulk file type:", bulk_type))

# Load individual files from the directory with filter suffix
# Use the specified bulk file type instead of always using "bulk"
bulk_path <- file.path(input_dir, paste0(dataset_prefix, "_", bulk_type, ".rda"))
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

# Determine the correct variable name for the bulk data
# This handles different naming conventions in different bulk files
if (exists("bulk")) {
  print("Found variable 'bulk' in the bulk file")
  bulk_data <- bulk
} else if (exists("bulk_random")) {
  print("Found variable 'bulk_random' in the bulk file")
  bulk_data <- bulk_random
} else if (exists("pseudobulk")) {
  print("Found variable 'pseudobulk' in the bulk file")
  bulk_data <- pseudobulk
} else {
  # Try to find any matrix in the loaded environment
  env_vars <- ls()
  matrix_vars <- env_vars[sapply(env_vars, function(v) is.matrix(get(v)))]
  
  if (length(matrix_vars) > 0) {
    bulk_data <- get(matrix_vars[1])
    print(paste("Using matrix variable:", matrix_vars[1]))
  } else {
    stop("Could not find bulk data matrix in the loaded file")
  }
}

print("CHECK: Data loaded successfully")

# Get common genes between bulk and single-cell data
common_genes <- intersect(rownames(bulk_data), rownames(singleCellExpr))
message(paste("Common genes between datasets:", length(common_genes)))
# Subset both matrices to common genes
bulk_data <- bulk_data[common_genes, , drop=FALSE]
singleCellExpr <- singleCellExpr[common_genes, , drop=FALSE]

# Get number of cell types from the single-cell data
n_cell_types <- length(unique(singleCellLabels))
print(paste("Number of cell types detected:", n_cell_types))

# Run CDSeq
print("Starting CDSeq deconvolution...")
cdseq.result <- CDSeq(
  bulk_data = bulk_data,
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

# Get the estimated proportions
prop_matrix <- t(cdseq.result$estProp)  # TRANSPOSED: Now samples in rows, cell types in columns

# Structure the results
deconvolutionResult <- list()
deconvolutionResult$CDSeq <- list(
  P = prop_matrix,
  S = cdseq.result$estGEP       # Genes in rows and cell types in columns
)

# Save results as RDA with filter suffix and bulk type indicated
results_filename <- file.path(output_dir, paste0("results_CDSeq_", sample_filter, "_", bulk_type, ".rda"))
save(deconvolutionResult, file=results_filename, compress=TRUE)

# Save proportions as CSV with filter suffix and bulk type indicated
proportions_df <- as.data.frame(deconvolutionResult$CDSeq$P)
proportions_df$Sample <- rownames(proportions_df)
csv_filename <- file.path(output_dir, paste0("results_CDSeq_", sample_filter, "_", bulk_type, "_proportions.csv"))
write.csv(proportions_df, file=csv_filename, row.names=FALSE)

print(paste("Results saved to:", results_filename))
print(paste("Proportions CSV saved to:", csv_filename))