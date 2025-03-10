#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
  stop("3 arguments must be supplied: R_LIBRARY_PATH INPUT_DIR OUTPUT_BASE_DIR", call.=FALSE)
}

####### Parameters
path_Rlibrary <- args[1]  # R library path
input_dir <- args[2]      # Directory containing the generated data files
output_base_dir <- args[3] # Base output directory for results

# Libraries
.libPaths(path_Rlibrary, FALSE)

# Load required packages
library(Biobase)  # For ExpressionSet creation
print("CHECK: Base libraries loaded")

# Install MuSiC if needed
if (!requireNamespace("MuSiC", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools", repos="https://cloud.r-project.org", lib=path_Rlibrary)
  }
  devtools::install_github("xuranw/MuSiC", lib=path_Rlibrary)
}
library(MuSiC)
print("CHECK: MuSiC library loaded")

# Extract dataset prefix from the input directory name
dataset_prefix <- basename(input_dir)

# Create dataset-specific output directory
output_dir <- file.path(output_base_dir, dataset_prefix)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
print(paste("Results will be saved to:", output_dir))

# Load individual files from the directory
bulk_path <- file.path(input_dir, paste0(dataset_prefix, "_bulk.rda"))
singleCellExpr_path <- file.path(input_dir, paste0(dataset_prefix, "_singleCellExpr.rda"))
singleCellLabels_path <- file.path(input_dir, paste0(dataset_prefix, "_singleCellLabels.rda"))

if (!file.exists(bulk_path)) {
  stop(paste("Bulk data file not found at:", bulk_path))
}
if (!file.exists(singleCellExpr_path)) {
  stop(paste("Single cell expression file not found at:", singleCellExpr_path))
}
if (!file.exists(singleCellLabels_path)) {
  stop(paste("Single cell labels file not found at:", singleCellLabels_path))
}

# Load the data
load(bulk_path)
load(singleCellExpr_path)
load(singleCellLabels_path)

print("CHECK: Data loaded successfully")
print(paste("Bulk dimensions:", paste(dim(bulk), collapse=" x ")))
print(paste("scRNA dimensions:", paste(dim(singleCellExpr), collapse=" x ")))

# Get common genes between bulk and single-cell data
common_genes <- intersect(rownames(bulk), rownames(singleCellExpr))
message(paste("Common genes between datasets:", length(common_genes)))

# Subset both matrices to common genes
bulk_filtered <- bulk[common_genes, , drop=FALSE]
sc_expr_filtered <- singleCellExpr[common_genes, , drop=FALSE]

print("Data dimensions after filtering to common genes:")
print(paste("Bulk:", paste(dim(bulk_filtered), collapse=" x ")))
print(paste("scRNA:", paste(dim(sc_expr_filtered), collapse=" x ")))

# Load optional single cell subject IDs if available
singleCellSubjects <- NULL
subjects_file <- file.path(input_dir, paste0(dataset_prefix, "_singleCellSubjects.rda"))
if (file.exists(subjects_file)) {
  load(subjects_file)
  print("Single cell subject IDs loaded successfully")
  print(paste("Number of subjects:", length(unique(singleCellSubjects))))
} else {
  # Create default subjects if not available (one subject for all cells)
  singleCellSubjects <- rep("subject1", length(singleCellLabels))
  print("Created default subject IDs - all cells assigned to one subject")
}

# Create single-cell metadata
sc_meta <- data.frame(
  SubjectName = singleCellSubjects,
  cellType = singleCellLabels,
  row.names = colnames(sc_expr_filtered)
)

# Create a SingleCellExperiment object for MuSiC
# First, create an ExpressionSet
sc_eset <- Biobase::ExpressionSet(
  assayData = sc_expr_filtered,
  phenoData = Biobase::AnnotatedDataFrame(sc_meta)
)

print("CHECK: ExpressionSet objects created")
print(paste("sc_eset dimensions:", paste(dim(Biobase::exprs(sc_eset)), collapse=" x ")))

# Run MuSiC deconvolution
print("Starting MuSiC deconvolution...")
music_result <- tryCatch({
  MuSiC::music_prop(
    bulk.mtx = bulk_filtered,    # Use the matrix directly
    sc.eset = sc_eset,           # Use the ExpressionSet for single cell
    clusters = "cellType",       # Column name for cell types
    samples = "SubjectName",     # Column name for subjects
    verbose = TRUE
  )
}, error = function(e) {
  print(paste("ERROR in MuSiC deconvolution:", e$message))
  return(NULL)
})

# Process and save results
if (!is.null(music_result)) {
  # Structure in the expected format
  deconvolutionResult <- list()
  deconvolutionResult$MuSiC <- list(
    P = music_result$Est.prop.weighted,  # Make sure this matches the format expected
    S = NULL             # Signature matrix not explicitly provided by MuSiC
  )
  
  # Save results
  music_results_filename <- file.path(output_dir, "results_MuSiC1.rda")
  save(deconvolutionResult, file=music_results_filename, compress=TRUE)
  
  # Save as CSV
  music_csv_filename <- file.path(output_dir, "results_MuSiC1_proportions.csv")
  music_props_df <- as.data.frame(music_result$Est.prop.weighted)
  music_props_df$Sample <- rownames(music_props_df)
  write.csv(music_props_df, file=music_csv_filename, row.names=FALSE)
  
  print(paste("MuSiC results saved to:", music_results_filename, "and", music_csv_filename))
} else {
  print("MuSiC results not saved due to error")
}

print("Deconvolution process completed")