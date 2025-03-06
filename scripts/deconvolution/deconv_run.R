#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
  stop(paste("5 arguments must be supplied instead of", length(args)), call. = FALSE)
}

####### Parameter of script (ORDER IS IMPORTANT)
path_Rlibrary <- args[1] #IMPORTANT
dataset_prefix <- args[2] # Dataset prefix/subfolder name
output_base_dir <- args[3] # Base output directory (deconv_results)
deconv_methods <- args[4] # Deconvolution method(s) - can be comma-separated list of methods
sparse_conversion <- as.logical(args[5])

# Libraries
.libPaths(path_Rlibrary, FALSE) #IMPORTANT
library(DeconBenchmark)
library(Matrix)
print("CHECK: Libraries loaded")

# Create dataset-specific output directory
output_dir <- file.path(output_base_dir, dataset_prefix)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
print(paste("Results will be saved to:", output_dir))

# Determine base path for loading files
base_data_dir <- file.path("/work/gr-fe/lorthiois/DeconBenchmark/generated_data", dataset_prefix)
print(paste("Loading data from:", base_data_dir))

# Load individual files
bulk_path <- file.path(base_data_dir, paste0(dataset_prefix, "_bulk.rda"))
singleCellExpr_path <- file.path(base_data_dir, paste0(dataset_prefix, "_singleCellExpr.rda"))
singleCellLabels_path <- file.path(base_data_dir, paste0(dataset_prefix, "_singleCellLabels.rda"))

# Check if files exist
if (!file.exists(bulk_path)) {
  stop(paste("Bulk data file not found:", bulk_path))
}
if (!file.exists(singleCellExpr_path)) {
  stop(paste("Single-cell expression file not found:", singleCellExpr_path))
}
if (!file.exists(singleCellLabels_path)) {
  stop(paste("Single-cell labels file not found:", singleCellLabels_path))
}

# Load data from individual files
load(bulk_path)
load(singleCellExpr_path)
load(singleCellLabels_path)

print("CHECK: Data loaded successfully")

# Check data dimensions for verification
print(paste("Bulk data dimensions:", paste(dim(bulk), collapse=" x ")))
print(paste("singleCellExpr dimensions:", paste(dim(singleCellExpr), collapse=" x ")))
print(paste("singleCellLabels length:", length(singleCellLabels)))

# Load additional optional files if they exist
markers <- NULL
markers_file <- file.path(base_data_dir, paste0(dataset_prefix, "_markers.rda"))
if (file.exists(markers_file)) {
  load(markers_file)
  print("Marker genes loaded successfully")
}

significantGenes <- NULL
sig_genes_file <- file.path(base_data_dir, paste0(dataset_prefix, "_significant_genes.rda"))
if (file.exists(sig_genes_file)) {
  load(sig_genes_file)
  print("Significant genes loaded successfully")
}

cellTypeExpr <- NULL
celltype_expr_file <- file.path(base_data_dir, paste0(dataset_prefix, "_celltype_expression.rda"))
if (file.exists(celltype_expr_file)) {
  load(celltype_expr_file)
  print("Cell type expression loaded successfully")
}

singleCellSubjects <- NULL
subjects_file <- file.path(base_data_dir, paste0(dataset_prefix, "_single_cell_subjects.rda"))
if (file.exists(subjects_file)) {
  load(subjects_file)
  print("Single cell subject IDs loaded successfully")
} else {
  # Create default subjects if not available
  singleCellSubjects <- rep("subject1", length(singleCellLabels))
  print("Created default subject IDs")
}

# Split method names
method_list <- unlist(strsplit(deconv_methods, ","))
method_list <- trimws(method_list)
print(paste("Methods to run:", paste(method_list, collapse=", ")))

# Generate references
reference <- generateReference(singleCellExpr, singleCellLabels, type="signature")
signature <- reference$signature

# # Generate missing inputs if needed
# if (is.null(cellTypeExpr) && any(grep("cellTypeExpr", unlist(getMethodsInputs(method_list))))) {
#   print("Auto-generating cell type expression matrix...")
#   # Code to generate cellTypeExpr
#   cell_types <- unique(singleCellLabels)
#   cellTypeExpr <- t(sapply(rownames(singleCellExpr), function(gene) {
#     sapply(cell_types, function(ct) {
#       mean(singleCellExpr[gene, singleCellLabels == ct])
#     })
#   }))
#   colnames(cellTypeExpr) <- cell_types
# }

# if (is.null(markers) && any(grep("markers", unlist(getMethodsInputs(method_list))))) {
#   print("Auto-generating marker genes...")
#   # Create a simple version of markers based on differential expression
#   markers <- list()
#   cell_types <- unique(singleCellLabels)
  
#   for (ct in cell_types) {
#     # Basic approach to find markers
#     cells_this_type <- singleCellLabels == ct
#     mean_expr_this <- rowMeans(singleCellExpr[, cells_this_type, drop=FALSE])
#     mean_expr_other <- rowMeans(singleCellExpr[, !cells_this_type, drop=FALSE])
    
#     # Calculate fold change
#     fold_change <- mean_expr_this / (mean_expr_other + 0.1)
    
#     # Get top genes by fold change
#     top_genes <- names(sort(fold_change, decreasing=TRUE)[1:min(20, length(fold_change))])
#     markers[[ct]] <- top_genes
#   }
# }

# if (is.null(significantGenes) && any(grep("significantGenes", unlist(getMethodsInputs(method_list))))) {
#   print("Using all genes as significant genes...")
#   significantGenes <- rownames(singleCellExpr)
# }

### Run deconvolution
print(paste("Starting deconvolution with", length(method_list), "methods..."))
deconvolutionResult <- runDeconvolution(
  methods = method_list,
  bulk = bulk,
  singleCellExpr = singleCellExpr,
  singleCellLabels = singleCellLabels,
  signature = signature,
  markers = markers,
  cellTypeExpr = cellTypeExpr,
  nCellTypes = length(unique(singleCellLabels)),
  singleCellSubjects = singleCellSubjects,
  isMethylation = FALSE,
  matlabLicenseFile = "303238", # MATLAB license
  containerEngine = "singularity"
)
print("CHECK: Deconvolution completed for all methods")

# Save results as before
# Create unified output with all methods
combined_method_name <- paste(method_list, collapse="_")
results_filename <- file.path(output_dir, paste0("results_", combined_method_name, ".rda"))
save(deconvolutionResult, file=results_filename)
print(paste("Results for all methods saved to:", results_filename))

# Also save individual method results for backward compatibility
for (method in method_list) {
  # Create a new list with the same structure as deconvolutionResult
  # but containing only the current method
  single_method_result <- list()
  single_method_result[[method]] <- deconvolutionResult[[method]]
  
  # Use the same variable name as the combined results
  temp_deconvolutionResult <- single_method_result
  
  single_results_filename <- file.path(output_dir, paste0("results_", method, ".rda"))
  save(temp_deconvolutionResult, file=single_results_filename, compress=TRUE)
  print(paste("Results for", method, "saved to:", single_results_filename))
}