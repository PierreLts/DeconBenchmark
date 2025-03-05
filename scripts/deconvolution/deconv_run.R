#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
  stop(paste("5 arguments must be supplied instead of", length(args)), call. = FALSE)
}

####### Parameter of script
path_Rlibrary <- args[1] #IMPORTANT
input_data <- args[2]
output_data <- args[3]
deconv_methods <- args[4] # can be comma-separated list of methods
sparse_conversion <- as.logical(args[5])

# Libraries
.libPaths(path_Rlibrary, FALSE) #IMPORTANT
library(DeconBenchmark)
library(Matrix)
print("CHECK: Libraries loaded")

# Load input data
loaded_objects <- load(input_data)
data_object_name <- loaded_objects[1]
data_object <- get(data_object_name)
bulk <- data_object$bulk
singleCellExpr <- data_object$singleCellExpr
singleCellLabels <- data_object$singleCellLabels
print("CHECK: Data extracted successfully")

# Load additional required inputs
markers <- NULL
markers_file <- file.path(dirname(input_data), "markers.rda")
if (file.exists(markers_file)) {
  load(markers_file)
  print("Marker genes loaded successfully")
}

significantGenes <- NULL
sig_genes_file <- file.path(dirname(input_data), "significant_genes.rda")
if (file.exists(sig_genes_file)) {
  load(sig_genes_file)
  print("Significant genes loaded successfully")
}

cellTypeExpr <- NULL
celltype_expr_file <- file.path(dirname(input_data), "celltype_expression.rda")
if (file.exists(celltype_expr_file)) {
  load(celltype_expr_file)
  print("Cell type expression loaded successfully")
}

singleCellSubjects <- NULL
subjects_file <- file.path(dirname(input_data), "single_cell_subjects.rda")
if (file.exists(subjects_file)) {
  load(subjects_file)
  print("Single cell subject IDs loaded successfully")
} else {
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

# Generate missing inputs if needed
if (is.null(cellTypeExpr)) {
  print("Generating cell type expression matrix...")
  # Code to generate cellTypeExpr from singleCellExpr and singleCellLabels
  cell_types <- unique(singleCellLabels)
  cellTypeExpr <- t(sapply(rownames(singleCellExpr), function(gene) {
    sapply(cell_types, function(ct) {
      mean(singleCellExpr[gene, singleCellLabels == ct])
    })
  }))
  colnames(cellTypeExpr) <- cell_types
}

if (is.null(markers)) {
  print("Generating marker genes...")
  markers <- list()
  cell_types <- unique(singleCellLabels)
  
  for (ct in cell_types) {
    # Get fold change between this cell type and others
    cells_this_type <- which(singleCellLabels == ct)
    cells_other_types <- which(singleCellLabels != ct)
    
    # Calculate mean expression
    mean_expr_this_type <- rowMeans(singleCellExpr[, cells_this_type, drop=FALSE])
    mean_expr_other_types <- rowMeans(singleCellExpr[, cells_other_types, drop=FALSE])
    
    # Calculate simple fold change
    fold_change <- mean_expr_this_type / (mean_expr_other_types + 0.1)
    
    # Get top genes 
    top_genes <- names(sort(fold_change, decreasing=TRUE)[1:20])
    markers[[ct]] <- top_genes
  }
}

if (is.null(significantGenes)) {
  print("Using all genes as significant genes...")
  significantGenes <- rownames(singleCellExpr)
}

### Run deconvolution
print(paste("Starting deconvolution with", length(method_list), "methods..."))
deconvolutionResult <- runDeconvolution(
  methods = method_list,
  bulk = bulk,
  singleCellExpr = singleCellExpr,
  singleCellLabels = singleCellLabels,
  signature = signature,
  markers = markers,
  #significantGenes = significantGenes,
  cellTypeExpr = cellTypeExpr,
  nCellTypes = length(unique(singleCellLabels)),
  singleCellSubjects = singleCellSubjects,
  isMethylation = FALSE,
  matlabLicenseFile = "303238", # MATLAB license
  containerEngine = "singularity"
)
print("CHECK: Deconvolution completed for all methods")

# Save results as before
# [Existing code for saving results]

# Print results preview
for (method in method_list) {
  if (!is.null(deconvolutionResult[[method]])) {
    print(paste("=== Results preview for", method, "==="))
    
    # Print available components
    components <- names(deconvolutionResult[[method]])
    print(paste("Available result components:", paste(components, collapse=", ")))
    
    # Preview proportions (P matrix)
    if (!is.null(deconvolutionResult[[method]]$P)) {
      proportion <- deconvolutionResult[[method]]$P
      print("Cell type proportions (P) preview:")
      print(head(proportion, 3))
    } else {
      print("No proportion matrix found")
    }

    # Preview signature matrix if available
    if (!is.null(deconvolutionResult[[method]]$S)) {
      signature <- deconvolutionResult[[method]]$S
      print("Signature matrix (S) preview:")
      # Handle different dimensions appropriately
      if (is.matrix(signature)) {
        print(head(signature[, 1:min(3, ncol(signature))], 3))
      } else if (is.vector(signature)) {
        print(head(signature, 6))
      } else {
        print("Signature matrix has unexpected format")
      }
    }
    
    # Preview other common outputs
    other_components <- setdiff(components, c("P", "S"))
    for (comp in other_components) {
      result_comp <- deconvolutionResult[[method]][[comp]]
      
      # Only print matrices or vectors with proper preview
      if (is.matrix(result_comp) || is.vector(result_comp)) {
        print(paste(comp, "preview:"))
        if (is.matrix(result_comp)) {
          print(head(result_comp[, 1:min(3, ncol(result_comp))], 3))
        } else if (length(result_comp) > 6) {
          print(result_comp[1:6])
        } else {
          print(result_comp)
        }
      } else {
        print(paste(comp, "is available but not previewed (non-tabular format)"))
      }
    }
    
    print("") # Empty line for separation
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