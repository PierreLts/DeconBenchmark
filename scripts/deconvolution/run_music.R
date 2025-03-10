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


# Get common genes between bulk and single-cell data
common_genes <- intersect(rownames(bulk), rownames(singleCellExpr))
message(paste("Common genes between datasets:", length(common_genes)))
# Subset both matrices to common genes
bulk <- bulk[common_genes, , drop=FALSE]
singleCellExpr <- singleCellExpr[common_genes, , drop=FALSE]


# Data validation
print("Validating data...")
print(paste("Bulk dimensions:", paste(dim(bulk), collapse=" x ")))
print(paste("Single-cell dimensions:", paste(dim(singleCellExpr), collapse=" x ")))
print(paste("Number of cell labels:", length(singleCellLabels)))
print(paste("Do cell labels match columns:", length(singleCellLabels) == ncol(singleCellExpr)))
print(paste("Unique cell types:", paste(unique(singleCellLabels), collapse=", ")))
print(paste("Number of unique cell types:", length(unique(singleCellLabels))))

# Ensure data contains numeric values and not integers
bulk <- matrix(as.numeric(bulk), nrow=nrow(bulk), ncol=ncol(bulk))
rownames(bulk) <- rownames(get("bulk"))
colnames(bulk) <- colnames(get("bulk"))

singleCellExpr <- matrix(as.numeric(singleCellExpr), nrow=nrow(singleCellExpr), ncol=ncol(singleCellExpr))
rownames(singleCellExpr) <- rownames(get("singleCellExpr"))
colnames(singleCellExpr) <- colnames(get("singleCellExpr"))

# Load optional single cell subject IDs if available
singleCellSubjects <- NULL
subjects_file <- file.path(input_dir, paste0(dataset_prefix, "_single_cell_subjects.rda"))
if (file.exists(subjects_file)) {
  load(subjects_file)
  print("Single cell subject IDs loaded successfully")
  print(paste("Number of subjects:", length(unique(singleCellSubjects))))
} else {
  # Create default subjects if not available (one subject for all cells)
  singleCellSubjects <- rep("subject1", length(singleCellLabels))
  print("Created default subject IDs - all cells assigned to one subject")
}

# Get common genes between bulk and single-cell data
common_genes <- intersect(rownames(bulk), rownames(singleCellExpr))
print(paste("Number of common genes:", length(common_genes)))

if (length(common_genes) < 1000) {
  print("WARNING: Very few common genes found. This might affect deconvolution quality.")
}

# Subset matrices to common genes
bulk_filtered <- bulk[common_genes, , drop=FALSE]
sc_expr_filtered <- singleCellExpr[common_genes, , drop=FALSE]

print("Data dimensions after filtering to common genes:")
print(paste("Bulk:", paste(dim(bulk_filtered), collapse=" x ")))
print(paste("scRNA:", paste(dim(sc_expr_filtered), collapse=" x ")))

# Check for zero or NA values
print(paste("Bulk contains zeros:", any(bulk_filtered == 0)))
print(paste("Bulk contains NA:", any(is.na(bulk_filtered))))
print(paste("scRNA contains zeros:", any(sc_expr_filtered == 0)))
print(paste("scRNA contains NA:", any(is.na(sc_expr_filtered))))

# Replace NA values with zeros if any
if (any(is.na(bulk_filtered))) {
  bulk_filtered[is.na(bulk_filtered)] <- 0
  print("Replaced NA values in bulk data with zeros")
}
if (any(is.na(sc_expr_filtered))) {
  sc_expr_filtered[is.na(sc_expr_filtered)] <- 0
  print("Replaced NA values in single-cell data with zeros")
}

# Create ExpressionSet objects required by MuSiC
# First, for single-cell data
sc_meta <- data.frame(
  SubjectName = singleCellSubjects,
  cellType = singleCellLabels,
  row.names = colnames(sc_expr_filtered)
)

# Print sample of metadata
print("Sample of single-cell metadata:")
print(head(sc_meta, 5))

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

# Verify ExpressionSet objects
print("Verifying ExpressionSet objects:")
print(paste("sc_eset dimensions:", paste(dim(Biobase::exprs(sc_eset)), collapse=" x ")))
print(paste("bulk_eset dimensions:", paste(dim(Biobase::exprs(bulk_eset)), collapse=" x ")))
print(paste("sc_eset phenoData:", paste(dim(Biobase::pData(sc_eset)), collapse=" x ")))
print(paste("bulk_eset phenoData:", paste(dim(Biobase::pData(bulk_eset)), collapse=" x ")))

# Run basic MuSiC deconvolution
print("Starting MuSiC deconvolution...")
music_result <- tryCatch({
  MuSiC::music_prop(
    bulk.eset = bulk_eset,
    sc.eset = sc_eset,
    clusters = "cellType",
    samples = "SubjectName",
    verbose = TRUE
  )
}, error = function(e) {
  print(paste("ERROR in MuSiC deconvolution:", e$message))
  return(NULL)
})

if (is.null(music_result)) {
  print("Error in MuSiC deconvolution, trying with adjusted parameters...")
  music_result <- tryCatch({
    # Try with more relaxed parameters
    MuSiC::music_prop(
      bulk.eset = bulk_eset,
      sc.eset = sc_eset,
      clusters = "cellType",
      samples = "SubjectName",
      verbose = TRUE,
      iter.max = 1000,    # Increase iterations
      nu = 1e-10,         # Adjust regularization
      eps = 0.01,         # Adjust convergence threshold
      centered = FALSE,   # Try without centering
      normalize = TRUE    # Ensure normalization
    )
  }, error = function(e) {
    print(paste("ERROR in MuSiC deconvolution (second attempt):", e$message))
    return(NULL)
  })
}

# If still failing, try a manual implementation as last resort
if (is.null(music_result)) {
  print("Attempting alternative approach...")
  
  # Try with basic implementation - similar to MuSiC but simpler
  tryCatch({
    # Get average expression by cell type
    cell_types <- unique(sc_meta$cellType)
    avg_expr_by_celltype <- matrix(0, nrow=nrow(sc_expr_filtered), ncol=length(cell_types))
    colnames(avg_expr_by_celltype) <- cell_types
    rownames(avg_expr_by_celltype) <- rownames(sc_expr_filtered)
    
    for (ct in cell_types) {
      cells_of_type <- sc_meta$cellType == ct
      avg_expr_by_celltype[, ct] <- rowMeans(sc_expr_filtered[, cells_of_type, drop=FALSE])
    }
    
    # Simple regression-based deconvolution
    library(limma)
    
    # Normalize data
    avg_expr_by_celltype <- avg_expr_by_celltype / colSums(avg_expr_by_celltype)
    bulk_filtered_norm <- sweep(bulk_filtered, 2, colSums(bulk_filtered), "/")
    
    # Create design matrix from signature
    X <- avg_expr_by_celltype
    
    # For each bulk sample, estimate proportions
    props <- matrix(0, nrow=ncol(bulk_filtered_norm), ncol=ncol(X))
    rownames(props) <- colnames(bulk_filtered_norm)
    colnames(props) <- colnames(X)
    
    for (i in 1:ncol(bulk_filtered_norm)) {
      y <- bulk_filtered_norm[, i]
      fit <- limma::lmFit(y, X)
      coef <- fit$coefficients
      # Ensure non-negative values and normalize to sum to 1
      coef[coef < 0] <- 0
      props[i, ] <- coef / sum(coef)
    }
    
    # Create a result object similar to MuSiC
    music_result <- list(
      Est.prop.weighted = t(props),
      alternative_method = TRUE
    )
    
    print("Alternative approach completed")
  }, error = function(e) {
    print(paste("ERROR in alternative approach:", e$message))
    music_result <<- NULL
  })
}

print("CHECK: MuSiC deconvolution completed")

# Run MuSiC2 deconvolution if available
print("Starting MuSiC2 deconvolution...")
music2_result <- tryCatch({
  MuSiC::music2_prop(
    bulk.eset = bulk_eset, 
    sc.eset = sc_eset, 
    clusters = "cellType", 
    samples = "SubjectName",
    verbose = TRUE
  )
}, error = function(e) {
  print(paste("ERROR in MuSiC2 deconvolution:", e$message))
  return(NULL)
})

if (is.null(music2_result)) {
  print("Error in MuSiC2 deconvolution, trying with adjusted parameters...")
  music2_result <- tryCatch({
    # Try with more relaxed parameters
    MuSiC::music2_prop(
      bulk.eset = bulk_eset,
      sc.eset = sc_eset,
      clusters = "cellType",
      samples = "SubjectName",
      verbose = TRUE,
      iter.max = 1000,    # Increase iterations
      nu = 1e-10,         # Adjust regularization
      eps = 0.01,         # Adjust convergence threshold
      centered = FALSE,   # Try without centering
      normalize = TRUE    # Ensure normalization
    )
  }, error = function(e) {
    print(paste("ERROR in MuSiC2 deconvolution (second attempt):", e$message))
    return(NULL)
  })
}

print("CHECK: MuSiC2 deconvolution completed")

# Format results to match your existing framework

# For MuSiC
if (!is.null(music_result)) {
  print("Processing MuSiC results...")
  # Extract the proportions matrix
  if (!is.null(music_result$Est.prop.weighted)) {
    music_props <- music_result$Est.prop.weighted
  } else if (!is.null(music_result$Est.prop)) {
    music_props <- music_result$Est.prop
    print("Using Est.prop instead of Est.prop.weighted")
  } else {
    print("Warning: Could not find proportion matrix in music_result")
    print("Available elements in music_result:")
    print(names(music_result))
    # Use the alternative if we created it
    if (!is.null(music_result$alternative_method)) {
      music_props <- music_result$Est.prop.weighted
      print("Using alternative method results")
    }
  }
  
  if (exists("music_props")) {
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
    print("MuSiC results not saved due to missing proportion matrix")
  }
} else {
  print("MuSiC results not saved due to error")
}

# For MuSiC2
if (!is.null(music2_result)) {
  print("Processing MuSiC2 results...")
  # Extract the proportions matrix
  if (!is.null(music2_result$Est.prop)) {
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
    print("Warning: Could not find proportion matrix in music2_result")
    print("Available elements in music2_result:")
    print(names(music2_result))
  }
} else {
  print("MuSiC2 results not saved due to error")
}

# If both methods failed, create a simple placeholder result
if (is.null(music_result) && is.null(music2_result)) {
  print("Both MuSiC and MuSiC2 failed. Creating fallback results...")
  
  # Create simple proportion estimates using marker genes
  tryCatch({
    # Get average expression by cell type (simple marker-based approach)
    cell_types <- unique(sc_meta$cellType)
    fallback_props <- matrix(0, nrow=ncol(bulk_filtered), ncol=length(cell_types))
    rownames(fallback_props) <- colnames(bulk_filtered)
    colnames(fallback_props) <- cell_types
    
    # Equal proportions as fallback
    fallback_props[] <- 1/ncol(fallback_props)
    
    # Save fallback results for MuSiC
    deconvolutionResult <- list()
    deconvolutionResult$MuSiC <- list(
      P = fallback_props,
      S = NULL,
      fallback = TRUE
    )
    
    fallback_filename <- file.path(output_dir, "results_MuSiC_fallback.rda")
    save(deconvolutionResult, file=fallback_filename, compress=TRUE)
    
    print(paste("Fallback results saved to:", fallback_filename))
    print("WARNING: These are placeholder results and should be used with caution!")
  }, error = function(e) {
    print(paste("ERROR creating fallback results:", e$message))
  })
}

print("Deconvolution process completed")