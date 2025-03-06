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

# Install required dependencies
print("Installing required dependencies...")

# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos="https://cloud.r-project.org", lib=path_Rlibrary)
}

# Install TOAST using BiocManager
if (!requireNamespace("TOAST", quietly = TRUE)) {
  BiocManager::install("TOAST", lib=path_Rlibrary)
}

# Install required CRAN packages
required_packages <- c("Matrix", "quadprog", "e1071", "devtools", "Biobase")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos="https://cloud.r-project.org", lib=path_Rlibrary)
  }
}

# Try to install MuSiC
print("Installing MuSiC package...")
if (!requireNamespace("MuSiC", quietly = TRUE)) {
  try_github <- try({
    devtools::install_github("xuranw/MuSiC", lib=path_Rlibrary)
  })
  
  # If GitHub install fails, try direct download of tarball
  if (inherits(try_github, "try-error")) {
    print("GitHub installation failed, trying direct download...")
    temp_dir <- tempdir()
    tarball_path <- file.path(temp_dir, "MuSiC.tar.gz")
    download.file("https://github.com/xuranw/MuSiC/archive/refs/tags/v1.0.0.tar.gz", 
                  destfile = tarball_path)
    install.packages(tarball_path, repos = NULL, type = "source", lib=path_Rlibrary)
  }
}

# Check if MuSiC was installed successfully
if (!requireNamespace("MuSiC", quietly = TRUE)) {
  stop("Failed to install MuSiC package. Please install it manually.")
}

# Load required packages
library(MuSiC)
library(Biobase)  # For ExpressionSet creation
print("CHECK: Libraries loaded")

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

# Load optional single cell subject IDs if available
singleCellSubjects <- NULL
subjects_file <- file.path(input_dir, paste0(dataset_prefix, "_single_cell_subjects.rda"))
if (file.exists(subjects_file)) {
  load(subjects_file)
  print("Single cell subject IDs loaded successfully")
} else {
  # Create default subjects if not available (one subject for all cells)
  singleCellSubjects <- rep("subject1", length(singleCellLabels))
  print("Created default subject IDs")
}

# Get common genes between bulk and single-cell data
common_genes <- intersect(rownames(bulk), rownames(singleCellExpr))
print(paste("Number of common genes:", length(common_genes)))

# Subset matrices to common genes
bulk_filtered <- bulk[common_genes, , drop=FALSE]
sc_expr_filtered <- singleCellExpr[common_genes, , drop=FALSE]

print("Data dimensions after filtering to common genes:")
print(paste("Bulk:", paste(dim(bulk_filtered), collapse=" x ")))
print(paste("scRNA:", paste(dim(sc_expr_filtered), collapse=" x ")))

# Create ExpressionSet objects required by MuSiC
# First, for single-cell data
sc_meta <- data.frame(
  SubjectName = singleCellSubjects,
  cellType = singleCellLabels,
  row.names = colnames(sc_expr_filtered)
)

# Create single-cell ExpressionSet
sc_eset <- Biobase::ExpressionSet(
  assayData = sc_expr_filtered,
  phenoData = Biobase::AnnotatedDataFrame(sc_meta)
)

# Create bulk ExpressionSet
bulk_meta <- data.frame(
  row.names = colnames(bulk_filtered)
)
bulk_eset <- Biobase::ExpressionSet(
  assayData = bulk_filtered,
  phenoData = Biobase::AnnotatedDataFrame(bulk_meta)
)

print("CHECK: ExpressionSet objects created")

# Run basic MuSiC deconvolution
print("Starting MuSiC deconvolution...")
music_result <- try({
  MuSiC::music_prop(
    bulk.eset = bulk_eset,
    sc.eset = sc_eset,
    clusters = "cellType",
    samples = "SubjectName",
    verbose = TRUE
  )
})

if (inherits(music_result, "try-error")) {
  print("Error in MuSiC deconvolution, trying with default parameters...")
  music_result <- try({
    MuSiC::music_prop(
      bulk.eset = bulk_eset,
      sc.eset = sc_eset,
      clusters = "cellType",
      samples = "SubjectName",
      verbose = TRUE,
      iter.max = 1000,  # Increase iterations
      nu = 1e-10,       # Adjust regularization
      eps = 0.01        # Adjust convergence threshold
    )
  })
}

print("CHECK: MuSiC deconvolution completed")

# Run MuSiC2 deconvolution if available
print("Starting MuSiC2 deconvolution...")
music2_result <- try({
  MuSiC::music2_prop(
    bulk.eset = bulk_eset, 
    sc.eset = sc_eset, 
    clusters = "cellType", 
    samples = "SubjectName",
    verbose = TRUE
  )
})

if (inherits(music2_result, "try-error")) {
  print("Error in MuSiC2 deconvolution, trying with default parameters...")
  music2_result <- try({
    MuSiC::music2_prop(
      bulk.eset = bulk_eset,
      sc.eset = sc_eset,
      clusters = "cellType",
      samples = "SubjectName",
      verbose = TRUE,
      iter.max = 1000,  # Increase iterations
      nu = 1e-10,       # Adjust regularization
      eps = 0.01        # Adjust convergence threshold
    )
  })
}

print("CHECK: MuSiC2 deconvolution completed")

# Format results to match your existing framework

# For MuSiC
if (!inherits(music_result, "try-error")) {
  # Extract the proportions matrix
  music_props <- music_result$Est.prop.weighted
  
  # Structure in the expected format
  deconvolutionResult <- list()
  deconvolutionResult$MuSiC <- list(
    P = t(music_props),  # Transpose to have samples in rows, cell types in columns
    S = NULL             # Signature matrix not explicitly provided by MuSiC
  )
  
  # Save results
  music_results_filename <- file.path(output_dir, "results_MuSiC.rda")
  save(deconvolutionResult, file=music_results_filename, compress=TRUE)
  
  # Save as CSV
  music_csv_filename <- file.path(output_dir, "results_MuSiC_proportions.csv")
  music_props_df <- as.data.frame(t(music_props))
  music_props_df$Sample <- rownames(music_props_df)
  write.csv(music_props_df, file=music_csv_filename, row.names=FALSE)
  
  print(paste("MuSiC results saved to:", music_results_filename, "and", music_csv_filename))
} else {
  print("MuSiC results not saved due to error")
}

# For MuSiC2
if (!inherits(music2_result, "try-error")) {
  # Extract the proportions matrix
  music2_props <- music2_result$Est.prop
  
  # Structure in the expected format
  deconvolutionResult <- list()
  deconvolutionResult$MuSiC2 <- list(
    P = t(music2_props),  # Transpose to have samples in rows, cell types in columns
    S = NULL              # Signature matrix not explicitly provided by MuSiC2
  )
  
  # Save results
  music2_results_filename <- file.path(output_dir, "results_MuSiC2.rda")
  save(deconvolutionResult, file=music2_results_filename, compress=TRUE)
  
  # Save as CSV
  music2_csv_filename <- file.path(output_dir, "results_MuSiC2_proportions.csv")
  music2_props_df <- as.data.frame(t(music2_props))
  music2_props_df$Sample <- rownames(music2_props_df)
  write.csv(music2_props_df, file=music2_csv_filename, row.names=FALSE)
  
  print(paste("MuSiC2 results saved to:", music2_results_filename, "and", music2_csv_filename))
} else {
  print("MuSiC2 results not saved due to error")
}

print("Deconvolution completed successfully")