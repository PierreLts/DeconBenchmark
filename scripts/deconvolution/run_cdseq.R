#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
  stop(paste("3 arguments must be supplied instead of", length(args)), call. = FALSE)
}

####### Parameter of script (ORDER IS IMPORTANT)
path_Rlibrary <- args[1] #IMPORTANT
input_data <- args[2]
output_data <- args[3]

# Libraries
.libPaths(path_Rlibrary, FALSE) #IMPORTANT
library(CDSeq)
print("CHECK: Libraries loaded")

# Load input data
loaded_objects <- load(input_data)
print("Loaded objects:")
print(loaded_objects)

# Get the 1st object
data_object_name <- loaded_objects[1]
# Access through the object name dynamically
data_object <- get(data_object_name)
bulk <- data_object$bulk
#singleCellExpr <- data_object$singleCellExpr
singleCellLabels <- data_object$singleCellLabels
print("CHECK: Data extracted successfully")

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

# # Cell type column naming
# if (ncol(deconvolutionResult$CDSeq$P) > 0) {
#   colnames(deconvolutionResult$CDSeq$P) <- paste0("CellType_", 1:ncol(deconvolutionResult$CDSeq$P))
# }

# Check if S needs to be transposed too
if (!is.null(deconvolutionResult$CDSeq$S) && nrow(deconvolutionResult$CDSeq$S) != nrow(bulk)) {
  print("Transposing S matrix to match expected format (genes in rows)")
  deconvolutionResult$CDSeq$S <- t(deconvolutionResult$CDSeq$S)
}

# Save results
results_filename <- file.path(output_data, paste0("results_CDSeq_Batch1.rda"))
save(deconvolutionResult, file=results_filename, compress=TRUE)
print(paste("Results saved to:", results_filename))