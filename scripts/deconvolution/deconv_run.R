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
# THE PIPELINE ONLY RUNS A METHOD AT ONCE
method_list <- unlist(strsplit(deconv_methods, ","))
method_list <- trimws(method_list)
method <- method_list[1]
print(paste("Methods to run:", paste(method_list, collapse=", ")))

# Generate references
reference <- generateReference(singleCellExpr, singleCellLabels, type="signature")
signature <- reference$signature



### Get common genes between bulk and single-cell data
common_genes <- intersect(rownames(bulk), rownames(singleCellExpr))
message(paste("Common genes between datasets:", length(common_genes)))

# Subset both matrices to common genes
bulk <- bulk[common_genes, , drop=FALSE]
singleCellExpr <- singleCellExpr[common_genes, , drop=FALSE]

# Verify dimensions match
message(paste("Filtered bulk matrix dimensions:", paste(dim(bulk), collapse=" x ")))
message(paste("Filtered scRNA matrix dimensions:", paste(dim(singleCellExpr), collapse=" x ")))
###

### Run deconvolution
#####
#####
# Applying fixes based on the specific method being run

print(paste("Applying pre-processing fixes for method:", method))

# For DeconPeaker: Thread count fix
if (method == "DeconPeaker") {
  print("Applying DeconPeaker thread count fix")
  # Set environment variable to limit threads to prevent "nthreads cannot be larger than NUMEXPR_MAX_THREADS" error
  Sys.setenv(NUMEXPR_MAX_THREADS="64")
}

# For DESeq2 and EMeth: Ensure gene dimensions match
if (method %in% c("DESeq2", "EMeth")) {
  print("Applying gene dimension fix for DESeq2/EMeth")
  if (!is.null(cellTypeExpr)) {
    common_genes <- intersect(rownames(bulk), rownames(cellTypeExpr))
    print(paste("Common genes between bulk and cellTypeExpr:", length(common_genes)))
    
    if (length(common_genes) > 0) {
      # Subset both matrices to common genes
      bulk <- bulk[common_genes, , drop=FALSE]
      cellTypeExpr <- cellTypeExpr[common_genes, , drop=FALSE]
      
      # Verify dimensions match
      print(paste("Fixed bulk matrix dimensions:", paste(dim(bulk), collapse=" x ")))
      print(paste("Fixed cellTypeExpr matrix dimensions:", paste(dim(cellTypeExpr), collapse=" x ")))
    } else {
      print("WARNING: No common genes found between bulk and cellTypeExpr!")
    }
  } else {
    print("WARNING: cellTypeExpr is NULL, cannot apply fix")
  }
}

# For DeCompress: Fix duplicate row names
if (method == "DeCompress") {
  print("Applying DeCompress duplicate rownames fix")
  
  # Check for duplicate rownames
  if (any(duplicated(rownames(bulk)))) {
    print(paste("Found", sum(duplicated(rownames(bulk))), "duplicate rownames in bulk data"))
    
    # Make rownames unique
    rownames(bulk) <- make.unique(rownames(bulk))
    
    # If cellTypeExpr exists, fix that too
    if (!is.null(cellTypeExpr) && any(duplicated(rownames(cellTypeExpr)))) {
      print(paste("Found", sum(duplicated(rownames(cellTypeExpr))), "duplicate rownames in cellTypeExpr"))
      rownames(cellTypeExpr) <- make.unique(rownames(cellTypeExpr))
    }
    
    print("Fixed duplicate rownames")
  } else {
    print("No duplicate rownames found, skipping fix")
  }
}

# For MIXTURE and EPIC: Ensure proper signature mapping
if (method %in% c("MIXTURE", "EPIC")) {
  print("Applying signature mapping fix for MIXTURE/EPIC")
  
  if (!is.null(signature)) {
    # Check if signature genes exist in bulk
    common_genes <- intersect(rownames(bulk), rownames(signature))
    print(paste("Common genes between bulk and signature:", length(common_genes)))
    
    if (length(common_genes) < 50) {
      print("WARNING: Few common genes found between bulk and signature!")
      
      # Try removing version numbers from Ensembl IDs if present
      if (any(grepl("\\.", rownames(bulk)))) {
        print("Attempting to remove Ensembl ID version numbers...")
        clean_rownames_bulk <- sub("\\..*", "", rownames(bulk))
        
        # Create temporary mapping from clean to original names
        temp_map <- setNames(rownames(bulk), clean_rownames_bulk)
        
        # Update bulk rownames temporarily for matching
        orig_rownames <- rownames(bulk)
        rownames(bulk) <- clean_rownames_bulk
        
        # Re-check common genes
        common_genes <- intersect(rownames(bulk), rownames(signature))
        print(paste("Common genes after Ensembl ID cleaning:", length(common_genes)))
        
        if (length(common_genes) >= 50) {
          # Keep only the common genes with cleaned IDs
          bulk <- bulk[common_genes, , drop=FALSE]
          signature <- signature[common_genes, , drop=FALSE]
          
          # Re-map bulk rownames back to original (but only for retained genes)
          remapped_rownames <- temp_map[rownames(bulk)]
          rownames(bulk) <- remapped_rownames
          
          print("Successfully fixed signature gene mapping")
        } else {
          # Revert to original rownames if fix didn't help
          rownames(bulk) <- orig_rownames
          print("Cleaning Ensembl IDs didn't improve gene matching")
        }
      }
    } else {
      # If enough common genes, just use those
      bulk <- bulk[common_genes, , drop=FALSE]
      signature <- signature[common_genes, , drop=FALSE]
      print("Signature gene mapping fixed")
    }
  } else {
    print("WARNING: signature is NULL, cannot apply fix")
  }
}

# Additional check for common genes between bulk and single-cell data
# Common gene matching for all methods - fundamental preprocessing step
common_genes <- intersect(rownames(bulk), rownames(singleCellExpr))
print(paste("Common genes between bulk and scRNA-seq data:", length(common_genes)))

# Only modify if we have some common genes but not all
if (length(common_genes) > 0 && 
    (length(common_genes) < nrow(bulk) || length(common_genes) < nrow(singleCellExpr))) {
  
  print("Applying gene matching between bulk and single-cell data")
  # Subset both matrices to common genes
  bulk <- bulk[common_genes, , drop=FALSE]
  singleCellExpr <- singleCellExpr[common_genes, , drop=FALSE]
  
  print(paste("Fixed bulk matrix dimensions:", paste(dim(bulk), collapse=" x ")))
  print(paste("Fixed scRNA matrix dimensions:", paste(dim(singleCellExpr), collapse=" x ")))
}

# Also match signature if it exists
if (!is.null(signature)) {
  common_genes_sig <- intersect(rownames(bulk), rownames(signature))
  if (length(common_genes_sig) > 0 && length(common_genes_sig) < nrow(signature)) {
    print("Applying gene matching between bulk and signature")
    bulk <- bulk[common_genes_sig, , drop=FALSE]
    signature <- signature[common_genes_sig, , drop=FALSE]
  }
}
#######
####### DECONVOLUTION #######
########################################################
# Now proceed with the deconvolution run
# Special handling for specific methods that need simpler parameters
if (method %in% c("TOAST", "RNA-Sieve", "PREDE", "BayICE", "AdRoit")) {
  print(paste("Applying simplified parameter set for method:", method))
  
  # Run with minimal parameter set
  deconvolutionResult <- tryCatch({
    runDeconvolution(
      methods = method,
      bulk = bulk,
      singleCellExpr = singleCellExpr,
      singleCellLabels = singleCellLabels,
      signature = signature,
      containerEngine = "singularity"
    )
  }, error = function(e) {
    # If this fails, print the error but return an empty result structure so the script can continue
    print(paste("ERROR running simplified method", method, ":", e$message))
    list(list(P = NULL, S = NULL, error = e$message))
  })
  
  print("CHECK: Simplified deconvolution completed")
} else {
  # Original code for other methods
  print("Using standard parameter set")
  deconvolutionResult <- runDeconvolution(
    methods = method,
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
}
######################################################



# Normalize deconvolution results (absolute values and sum to 1 per sample)
if (!is.null(deconvolutionResult[[method]]$P)) {
  prop_matrix <- deconvolutionResult[[method]]$P
  
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
  deconvolutionResult[[method]]$P <- prop_matrix
  print("Normalized cell type proportions (absolute values, sum to 1 per sample)")
}




# Save as RDA format
results_filename <- file.path(output_dir, paste0("results_", method, ".rda"))
save(deconvolutionResult, file=results_filename, compress=TRUE)
print(paste("Results saved to:", results_filename))

# Save as CSV format
# Extract the proportions matrix (P) and convert to dataframe
if (!is.null(deconvolutionResult[[method]]$P)) {
  proportions_df <- as.data.frame(deconvolutionResult[[method]]$P)
  
  # Add sample names as a column
  proportions_df$Sample <- rownames(proportions_df)
  
  # Save as CSV
  csv_filename <- file.path(output_dir, paste0("results_", method, "_proportions.csv"))
  write.csv(proportions_df, file=csv_filename, row.names=FALSE)
  print(paste("Proportions CSV saved to:", csv_filename))
}

# If signature matrix exists, save it too
if (!is.null(deconvolutionResult[[method]]$S)) {
  signature_df <- as.data.frame(deconvolutionResult[[method]]$S)
  
  # Add gene names as a column
  signature_df$Gene <- rownames(signature_df)
  
  # Save as CSV
  sig_csv_filename <- file.path(output_dir, paste0("results_", method, "_signature.csv"))
  write.csv(signature_df, file=sig_csv_filename, row.names=FALSE)
  print(paste("Signature CSV saved to:", sig_csv_filename))
}