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


# Get common genes between bulk and single-cell data
common_genes <- intersect(rownames(bulk), rownames(singleCellExpr))
message(paste("Common genes between datasets:", length(common_genes)))
# Subset both matrices to common genes
bulk <- bulk[common_genes, , drop=FALSE]
singleCellExpr <- singleCellExpr[common_genes, , drop=FALSE]

# Verify dimensions match
message(paste("Filtered bulk matrix dimensions:", paste(dim(bulk), collapse=" x ")))
message(paste("Filtered scRNA matrix dimensions:", paste(dim(singleCellExpr), collapse=" x ")))



# Split method names ()
method_list <- unlist(strsplit(deconv_methods, ","))
method_list <- trimws(method_list)
method <- method_list[1]
print(paste("Methods to run:", paste(method)))





# # Bulk normalization
# if (method %in% c("MuSiC","DWLS","AdRoit","spatialDWLS","scaden","DigitalDLSorter","RNA-Sieve","DecOT","MOMF","DeMixt","FARDEEP","ARIC","ImmuCellAI","EPIC","DeCompress","TICPE","Linseed","TOAST","BayCount","SMC","Deblender","MCPcounter","Bseq-SC","BayesPrism")) {
#   # Raw count
# }

# TPM
if (method %in% c("LinDeconSeq","BayICE","DESeq2")) {
  # TPM normalization
}

# Log transformed
if (method %in% c("dtangle","PREDE")) {
}

# Normalized transcriptional measurements
if (method %in% c("DeconRNASeq")) {
}


# other data
# Methylation data
if (method %in% c("BayesCCE","ReFACTor","EMeth")) {
}

# DeconPeaker needs ATAC-Seq & Box-Cox transformation
if (method %in% c("DeconPeaker")) {
}






# # For DeconPeaker: Thread count fix
# if (method == "DeconPeaker") {
#   print("Applying DeconPeaker thread count fix")
#   # Set environment variable to limit threads to prevent "nthreads cannot be larger than NUMEXPR_MAX_THREADS" error
#   Sys.setenv(NUMEXPR_MAX_THREADS="64")
# }



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
  singleCellExpr = singleCellExpr,
  singleCellLabels = singleCellLabels,
  signature = signature,
  markers = markers,
  cellTypeExpr = cellTypeExpr,
  nCellTypes = length(unique(singleCellLabels)),
  # singleCellSubjects = singleCellSubjects,
  isMethylation = FALSE,
  matlabLicenseFile = "303238", # MATLAB license
  containerEngine = "singularity"
)
print("CHECK: Deconvolution completed for all methods")



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