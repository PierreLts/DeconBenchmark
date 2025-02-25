#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
  stop(paste("5 arguments must be supplied instead of", length(args)), call. = FALSE)
}

####### Parameter of script (ORDER IS IMPORTANT)
path_Rlibrary <- args[1] #IMPORTANT
input_data <- args[2]
output_data <- args[3]
deconv_methods <- args[4] # Now can be comma-separated list of methods
sparse_conversion <- as.logical(args[5])

# Libraries
.libPaths(path_Rlibrary, FALSE) #IMPORTANT
library(DeconBenchmark)
library(Matrix)
print("CHECK: Libraries loaded")


# Set up custom image locations
singularity_images <- list(
  "bseqsc" = "/scratch/lorthiois/singularity_images/bseqsc.sif",
  "celldistinguisher" = "/scratch/lorthiois/singularity_images/celldistinguisher.sif",
  "cibersort" = "/scratch/lorthiois/singularity_images/cibersort.sif",
  "demixt" = "/scratch/lorthiois/singularity_images/demixt.sif",
  "methylresolver" = "/scratch/lorthiois/singularity_images/methylresolver.sif"
)

 
# Load input data
loaded_objects <- load(input_data)
print("Loaded objects:")
print(loaded_objects)

# Get the name of the loaded object (assuming it's the first/only one)
data_object_name <- loaded_objects[1]
# Access through the object name dynamically
data_object <- get(data_object_name)
bulk <- data_object$bulk
singleCellExpr <- data_object$singleCellExpr
singleCellLabels <- data_object$singleCellLabels
print("CHECK: Data extracted successfully")

if (sparse_conversion) {
  print("Converting to sparse matrices for memory efficiency...")
  # Convert to sparse format for memory efficiency
  sparse_singleCellExpr <- Matrix(singleCellExpr, sparse=TRUE)
  sparse_bulk <- Matrix(bulk, sparse=TRUE)
  
  # Clear originals to free memory
  rm(singleCellExpr, bulk)
  gc()
  
  # Convert back to dense only when passing to container
  singleCellExpr <- as.matrix(sparse_singleCellExpr)
  bulk <- as.matrix(sparse_bulk)
  
  # Clear sparse versions
  rm(sparse_singleCellExpr, sparse_bulk)
  gc()
  print("Memory optimization complete")
}

# Split the comma-separated method names
method_list <- unlist(strsplit(deconv_methods, ","))
method_list <- trimws(method_list) # Remove any whitespace
print(paste("Methods to run:", paste(method_list, collapse=", ")))

# Generate reference if needed for multiple methods
signature <- NULL
if (length(method_list) > 1) {
  # Check if any methods require a signature
  required_inputs <- getMethodsInputs(method_list, containerEngine = "singularity")
  needs_signature <- any(sapply(required_inputs, function(x) "signature" %in% x))
  
  if (needs_signature) {
    print("Generating reference signature for methods that require it...")
    reference <- generateReference(singleCellExpr, singleCellLabels, type="signature")
    signature <- reference$signature
  }
}

# Run deconvolution using Singularity with all specified methods
print(paste("Starting deconvolution with", length(method_list), "methods..."))
deconvolutionResult <- runDeconvolution(
  methods = method_list,
  bulk = bulk,
  singleCellExpr = singleCellExpr,
  singleCellLabels = singleCellLabels,
  signature = signature,
  containerEngine = "singularity"
)
print("CHECK: Deconvolution completed for all methods")

# Print results preview
for (method in method_list) {
  if (!is.null(deconvolutionResult[[method]])) {
    proportion <- deconvolutionResult[[method]]$P
    if (!is.null(proportion)) {
      print(paste("Preview of", method, "results:"))
      print(head(proportion, 3))
    } else {
      print(paste("No proportion matrix found for", method))
    }
  } else {
    print(paste("Method", method, "failed to produce results"))
  }
}

# Save results
# Extract the file name without the path
input_filename <- basename(input_data)
input_filename <- tools::file_path_sans_ext(input_filename) # Remove file extension

# Create unified output with all methods
results_filename <- file.path(output_data, paste0("results_", 
                                             paste(method_list, collapse="_"), 
                                             "_", input_filename, ".rda"))
save(deconvolutionResult, file=results_filename)
print(paste("Results for all methods saved to:", results_filename))

# Also save individual method results for backward compatibility
for (method in method_list) {
  # Create a new list with the same structure as deconvolutionResult
  # but containing only the current method
  single_method_result <- list()
  single_method_result[[method]] <- deconvolutionResult[[method]]
  
  # Use the same variable name as the combined results
  deconvolutionResult <- single_method_result
  
  single_results_filename <- file.path(output_data, paste0("results_", method, "_", input_filename, ".rda"))
  save(deconvolutionResult, file=single_results_filename, compress=TRUE)
  print(paste("Results for", method, "saved to:", single_results_filename))
}

# Restore the full results for the final output message
deconvolutionResult <- get("deconvolutionResult", inherits=FALSE)