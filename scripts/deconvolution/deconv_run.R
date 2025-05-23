#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 7) {
  stop(paste("7 arguments must be supplied instead of", length(args)), call. = FALSE)
}

####### Parameter of script (ORDER IS IMPORTANT)
path_Rlibrary <- args[1] #IMPORTANT
dataset_prefix <- args[2] # Dataset prefix/subfolder name
sample_filter <- args[3] # Sample filter: A, B, or AB
output_base_dir <- args[4] # Base output directory (deconv_results)
deconv_methods <- args[5] # Deconvolution method(s) - can be comma-separated list of methods
sparse_conversion <- as.logical(args[6])
bulk_type <- args[7] # Bulk file type (e.g., bulk, bulk_random, pseudobulk)


# Libraries
.libPaths(path_Rlibrary, FALSE) #IMPORTANT
library(DeconBenchmark)
library(Matrix)
print("CHECK: Libraries loaded")

# Create dataset-specific output directory
output_dir <- output_base_dir

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
print(paste("Results will be saved to:", output_dir))

# Determine base path for loading files
base_data_dir <- file.path("/work/gr-fe/lorthiois/DeconBenchmark/generated_data", dataset_prefix)
print(paste("Loading data from:", base_data_dir))
print(paste("Using sample filter:", sample_filter))

# Construct file paths with appropriate sample filter suffix
# Bulk data doesn't have a filter suffix
bulk_path <- file.path(base_data_dir, paste0(dataset_prefix, "_", bulk_type, ".rda"))
# Filtered data files
singleCellExpr_path <- file.path(base_data_dir, paste0(dataset_prefix, "_singleCellExpr_", sample_filter, ".rda"))
singleCellLabels_path <- file.path(base_data_dir, paste0(dataset_prefix, "_singleCellLabels_", sample_filter, ".rda"))

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

# Get common genes between bulk and single-cell data
common_genes <- intersect(rownames(bulk), rownames(singleCellExpr))
message(paste("Common genes between datasets:", length(common_genes)))
# Subset both matrices to common genes
bulk <- bulk[common_genes, , drop=FALSE]
singleCellExpr <- singleCellExpr[common_genes, , drop=FALSE]

# Verify dimensions match
message(paste("Filtered bulk matrix dimensions:", paste(dim(bulk), collapse=" x ")))
message(paste("Filtered scRNA matrix dimensions:", paste(dim(singleCellExpr), collapse=" x ")))
### Input data loaded and filtered ####

# Split method names ()
method_list <- unlist(strsplit(deconv_methods, ","))
method_list <- trimws(method_list)
method <- method_list[1]
print(paste("Methods to run:", paste(method)))

# Load optional single cell subject IDs if available
singleCellSubjects <- NULL
subjects_file <- file.path(base_data_dir, paste0(dataset_prefix, "_singleCellSubjects_", sample_filter, ".rda"))
if (file.exists(subjects_file)) {
  load(subjects_file)
  print("Single cell subject IDs loaded successfully")
} else {
  # Create default subjects if not available
  singleCellSubjects <- rep("subject1", length(singleCellLabels))
  print("Created default subject IDs")
}

### Normalization ###
# TPM
if (method %in% c("LinDeconSeq","BayICE","DESeq2")) {
  print('SHOULD APPLY TPM NORM')
  # TPM normalization
}

# Log transformed (If it means log2() ?)
if (method %in% c("dtangle","PREDE")) {
  print("Applying log transformation...")
  # log2(x+1) transformation
  pseudo_count <- 1
  bulk_log <- log2(bulk + pseudo_count)
  bulk <- bulk_log
  print("Bulk data log-transformed (log2(counts+1))")
}

# Generate references
reference <- generateReference(singleCellExpr, singleCellLabels, c("markers", "sigGenes", "signature", "cellTypeExpr"), 1)
markers <- reference$markers
sigGenes <- reference$sigGenes
signature <- reference$signature
cellTypeExpr <- reference$cellTypeExpr

####### DECONVOLUTION ##################################
print("Using standard parameter set")
deconvolutionResult <- runDeconvolution(
  methods = method,
  bulk = bulk,
  nCellTypes = length(unique(singleCellLabels)),
  markers = markers,
  isMethylation = FALSE,
  singleCellExpr = singleCellExpr,
  singleCellLabels = singleCellLabels,
  singleCellSubjects = singleCellSubjects, # BisqueRef
  cellTypeExpr = cellTypeExpr,
  sigGenes = sigGenes,
  signature = signature,
  seed = 1,
  matlabLicenseFile = "303238", # MATLAB license
  containerEngine = "singularity"
)
print("CHECK: Deconvolution completed for all methods")

# Normalize deconvolution results - all-or-nothing approach for entire dataset to convert percentage to 0-1
if (!is.null(deconvolutionResult[[method]]$P)) {
  prop_matrix <- deconvolutionResult[[method]]$P
  
  # Check if ALL samples have sums close to 100
  all_sums_near_100 <- TRUE
  for (i in 1:nrow(prop_matrix)) {
    row_sum <- sum(prop_matrix[i, ])
    if (abs(row_sum - 100) >= 1) {
      all_sums_near_100 <- FALSE
      print(paste("Sample", rownames(prop_matrix)[i], "has sum", row_sum, "- not applying percentage normalization"))
      break  # Stop checking once we find one that doesn't match
    }
  }
  
  # Apply normalization only if ALL samples sum to approximately 100
  if (all_sums_near_100) {
    # Normalize from percentage (0-100) to proportion (0-1)
    prop_matrix <- prop_matrix / 100
    print("All samples had sums close to 100, normalized from percentages (0-100) to proportions (0-1)")
  } else {
    print("Not all samples had sums close to 100, preserving original values")
  }
  
  # Update the proportions matrix in the result object
  deconvolutionResult[[method]]$P <- prop_matrix
}

# Save as RDA format with filter suffix
results_filename <- file.path(output_dir, paste0("results_", method, "_", sample_filter, ".rda"))
save(deconvolutionResult, file=results_filename, compress=TRUE)
print(paste("Results saved to:", results_filename))

# Save as CSV format with filter suffix
# Extract the proportions matrix (P) and convert to dataframe
if (!is.null(deconvolutionResult[[method]]$P)) {
  proportions_df <- as.data.frame(deconvolutionResult[[method]]$P)
  
  # Add sample names as a column
  proportions_df$Sample <- rownames(proportions_df)
  
  # Save as CSV
  csv_filename <- file.path(output_dir, paste0("results_", method, "_", sample_filter, "_proportions.csv"))
  write.csv(proportions_df, file=csv_filename, row.names=FALSE)
  print(paste("Proportions CSV saved to:", csv_filename))
}

# If signature matrix exists, save it too
if (!is.null(deconvolutionResult[[method]]$S)) {
  signature_df <- as.data.frame(deconvolutionResult[[method]]$S)
  
  # Add gene names as a column
  signature_df$Gene <- rownames(signature_df)
  
  # Save as CSV
  sig_csv_filename <- file.path(output_dir, paste0("results_", method, "_", sample_filter, "_signature.csv"))
  write.csv(signature_df, file=sig_csv_filename, row.names=FALSE)
  print(paste("Signature CSV saved to:", sig_csv_filename))
}