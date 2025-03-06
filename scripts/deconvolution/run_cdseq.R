#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
  stop(paste("3 arguments must be supplied instead of", length(args)), call. = FALSE)
}

####### Parameter of script (ORDER IS IMPORTANT)
path_Rlibrary <- args[1] #IMPORTANT
input_dir <- args[2]     # Directory containing the generated data files
output_base_dir <- args[3]   # Base output directory for results

# Libraries
.libPaths(path_Rlibrary, FALSE) #IMPORTANT
library(CDSeq)
print("CHECK: Libraries loaded")

# Extract dataset prefix from the input directory name
dataset_prefix <- basename(input_dir)

# Create dataset-specific output directory
output_dir <- file.path(output_base_dir, dataset_prefix)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
print(paste("Results will be saved to:", output_dir))

# Load individual files from the directory
bulk_path <- file.path(input_dir, paste0(dataset_prefix, "_bulk.rda"))
singleCellLabels_path <- file.path(input_dir, paste0(dataset_prefix, "_singleCellLabels.rda"))

if (!file.exists(bulk_path)) {
  stop(paste("Bulk data file not found at:", bulk_path))
}
if (!file.exists(singleCellLabels_path)) {
  stop(paste("Single cell labels file not found at:", singleCellLabels_path))
}

# Load the data
load(bulk_path)
load(singleCellLabels_path)

print("CHECK: Data loaded successfully")

# Get number of cell types from the single-cell data
n_cell_types <- length(unique(singleCellLabels))
print(paste("Number of cell types detected:", n_cell_types))

# Run CDSeq
print("Starting CDSeq deconvolution...")
cdseq_result <- CDSeq(
  bulk_data = bulk,
  cell_type_number = n_cell_types, #seq(max(2, n_cell_types-3), min(n_cell_types+3, 15), 3),  # Try a range around the known number
  beta = 0.5,                 
  alpha = 5,                           
  mcmc_iterations = 700,             # (700-2000)
  dilution_factor = 5,               # (2-10)
  gene_subset_size = 300,            # (200-500)
  block_number = 5,                  # (>5)
  cpu_number = 16                    # CPUs
)

print("CHECK: CDSeq deconvolution completed")

# Structure the results to match plotting format
# Transpose the matrices
deconvolutionResult <- list()
deconvolutionResult$CDSeq <- list(
  P = t(cdseq_result$estProp),  # TRANSPOSED: Now samples in rows, cell types in columns
  S = cdseq_result$estGEP       # Check if this has genes in rows and cell types in columns
)

# Save results
results_filename <- file.path(output_dir, "results_CDSeq.rda")
save(deconvolutionResult, file=results_filename, compress=TRUE)
print(paste("Results saved to:", results_filename))