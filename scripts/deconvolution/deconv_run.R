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

#########################################################
# GENERAL DATA PREPROCESSING - Applied to all methods
#########################################################

print("Applying general data preprocessing...")

# 1. Clean and standardize all gene IDs (removing version numbers)
clean_gene_ids <- function(mat) {
  if (is.null(mat)) return(NULL)
  # Remove Ensembl version numbers (e.g., ENSG00000123456.1 -> ENSG00000123456)
  new_rownames <- sub("\\.\\d+$", "", rownames(mat))
  # Check for duplicates after cleaning
  if (any(duplicated(new_rownames))) {
    print("WARNING: Duplicate gene IDs found after cleaning - making unique")
    new_rownames <- make.unique(new_rownames)
  }
  rownames(mat) <- new_rownames
  return(mat)
}

# Apply to all matrices
bulk <- clean_gene_ids(bulk)
singleCellExpr <- clean_gene_ids(singleCellExpr)
if (!is.null(cellTypeExpr)) cellTypeExpr <- clean_gene_ids(cellTypeExpr)

# 2. Handle NA values in all matrices
replace_na_with_zero <- function(mat) {
  if (is.null(mat)) return(NULL)
  na_count <- sum(is.na(mat))
  if (na_count > 0) {
    print(paste("Replacing", na_count, "NA values with zeros"))
    mat[is.na(mat)] <- 0
  }
  return(mat)
}

bulk <- replace_na_with_zero(bulk)
singleCellExpr <- replace_na_with_zero(singleCellExpr)
if (!is.null(cellTypeExpr)) cellTypeExpr <- replace_na_with_zero(cellTypeExpr)

# 3. Verify no zero-variance genes
remove_zero_var_genes <- function(mat) {
  if (is.null(mat)) return(NULL)
  var_by_gene <- apply(mat, 1, var)
  zero_var_idx <- which(var_by_gene == 0)
  if (length(zero_var_idx) > 0) {
    print(paste("Removing", length(zero_var_idx), "zero-variance genes"))
    mat <- mat[-zero_var_idx, , drop=FALSE]
  }
  return(mat)
}

bulk <- remove_zero_var_genes(bulk)
singleCellExpr <- remove_zero_var_genes(singleCellExpr)
if (!is.null(cellTypeExpr)) cellTypeExpr <- remove_zero_var_genes(cellTypeExpr)

# 4. Generate reference materials
all_cell_types <- unique(singleCellLabels)
n_cell_types <- length(all_cell_types)
print(paste("Number of cell types:", n_cell_types))
print(paste("Cell types:", paste(all_cell_types, collapse=", ")))

# Generate average expression per cell type (useful for many methods)
avg_expr_by_celltype <- matrix(0, nrow=nrow(singleCellExpr), ncol=length(all_cell_types))
rownames(avg_expr_by_celltype) <- rownames(singleCellExpr)
colnames(avg_expr_by_celltype) <- all_cell_types

for (ct in all_cell_types) {
  cells_of_type <- singleCellLabels == ct
  if (sum(cells_of_type) > 0) {
    avg_expr_by_celltype[, ct] <- rowMeans(singleCellExpr[, cells_of_type, drop=FALSE])
  }
}
print("Created average expression matrix by cell type")

# Generate robust marker genes if not available or insufficient
if (is.null(markers) || length(markers) < length(all_cell_types) || 
    any(sapply(markers, length) < 10)) {
  print("Generating robust marker genes")
  robust_markers <- list()
  
  # Calculate fold change vs. other cell types
  for (ct in all_cell_types) {
    expr_in_ct <- avg_expr_by_celltype[, ct]
    expr_in_others <- rowMeans(avg_expr_by_celltype[, colnames(avg_expr_by_celltype) != ct, drop=FALSE])
    
    # Simple log fold change
    log_fc <- log2((expr_in_ct + 1) / (expr_in_others + 1))
    specificity <- expr_in_ct / (expr_in_others + 1e-10)
    
    # Combined score weighting both fold change and specificity
    marker_score <- log_fc * specificity
    
    # Select top genes with highest score (minimum 20 per cell type)
    gene_rank <- rank(-marker_score)  # Negative to rank highest first
    min_markers <- min(50, sum(log_fc > 1 & specificity > 2))
    if (min_markers < 20) min_markers <- 20
    top_genes <- names(sort(marker_score, decreasing=TRUE)[1:min_markers])
    
    robust_markers[[ct]] <- top_genes
  }
  
  # Use both original markers and robust markers if original exists
  if (!is.null(markers)) {
    combined_markers <- markers
    for (ct in names(robust_markers)) {
      if (ct %in% names(markers)) {
        # Add new markers to existing ones
        combined_markers[[ct]] <- unique(c(markers[[ct]], robust_markers[[ct]]))
      } else {
        # New cell type
        combined_markers[[ct]] <- robust_markers[[ct]]
      }
    }
    markers <- combined_markers
  } else {
    markers <- robust_markers
  }
  
  print("Final marker counts per cell type:")
  marker_counts <- sapply(markers, length)
  print(marker_counts)
}

# Generate standard signature
reference <- generateReference(singleCellExpr, singleCellLabels, type="signature")
signature <- reference$signature

#########################################################
# METHOD-SPECIFIC PREPROCESSING
#########################################################

print(paste("Applying pre-processing fixes for method:", method))

# AdRoit - Fix mapper object issue
if (method == "AdRoit") {
  print("Enhanced AdRoit fix")
  # Create mapper global variable and save to multiple environments
  # This addresses issues of variable scope and visibility
  tryCatch({
    mapper <- data.frame(
      gene_id = rownames(singleCellExpr),
      gene_name = rownames(singleCellExpr),
      stringsAsFactors = FALSE
    )
    assign("mapper", mapper, envir = .GlobalEnv)
    
    # Try assigning to parent environment as well for visibility
    parent.env <- parent.frame()
    assign("mapper", mapper, envir = parent.env)
    
    # Copy to package namespace (as a desperate measure)
    tryCatch({
      ns <- asNamespace("DeconBenchmark")
      assign("mapper", mapper, envir = ns)
    }, error = function(e) {
      print("Note: Could not save mapper to DeconBenchmark namespace")
    })
    
    # Also try direct assignment in case the issue is with `assign`
    eval(parse(text = "mapper <- mapper"), envir = .GlobalEnv)
    
    print("AdRoit mapper created successfully")
    
    # Manually harmonize gene IDs for AdRoit
    common_genes <- intersect(rownames(bulk), rownames(singleCellExpr))
    bulk <- bulk[common_genes, , drop=FALSE]
    singleCellExpr <- singleCellExpr[common_genes, , drop=FALSE]
    
    # Update mapper to match the filtered genes
    mapper <- mapper[mapper$gene_id %in% common_genes, ]
    assign("mapper", mapper, envir = .GlobalEnv)
    
    # Write mapper to disk to ensure it survives container issues
    temp_mapper_file <- tempfile(pattern = "mapper_", fileext = ".rda")
    save(mapper, file = temp_mapper_file)
    # Force global variable
    .GlobalEnv$mapper <- mapper
  }, error = function(e) {
    print(paste("Warning: AdRoit mapper creation failed:", e$message))
  })
}

# BayICE - Fix matrix indexing
if (method == "BayICE") {
  print("Enhanced BayICE fix")
  # More robust gene matching and order verification
  common_genes <- intersect(rownames(bulk), rownames(singleCellExpr))
  print(paste("Common genes between bulk and scRNA:", length(common_genes)))
  
  # Check if number of common genes is reasonable
  if (length(common_genes) < 100) {
    stop("Too few common genes for BayICE to work properly")
  }
  
  # Sort for consistent ordering
  common_genes <- sort(common_genes)
  
  # Subset matrices
  bulk <- bulk[common_genes, , drop=FALSE]
  singleCellExpr <- singleCellExpr[common_genes, , drop=FALSE]
  
  # Verify all rownames match exactly - this is critical for BayICE
  if (!all(rownames(bulk) == rownames(singleCellExpr))) {
    print("WARNING: Forcing exact rowname match for BayICE")
    # Force exact match
    sc_expr_copy <- singleCellExpr
    rownames(sc_expr_copy) <- rownames(bulk)
    singleCellExpr <- sc_expr_copy
  }
  
  print("Matrices aligned and verified for BayICE")
  print(paste("Final dimensions - bulk:", paste(dim(bulk), collapse="x"), 
              "scRNA:", paste(dim(singleCellExpr), collapse="x")))
}

# PREDE - Fix feature selection
if (method == "PREDE") {
  print("Enhanced PREDE fix")
  # PREDE needs cellTypeExpr matrix with matching dimensions
  
  # If cellTypeExpr doesn't exist, create it from singleCellExpr
  if (is.null(cellTypeExpr)) {
    print("Creating cellTypeExpr from singleCellExpr for PREDE")
    cellTypeExpr <- avg_expr_by_celltype
  }
  
  # Align genes across all matrices
  common_genes <- Reduce(intersect, list(
    rownames(bulk), 
    rownames(singleCellExpr),
    rownames(cellTypeExpr)
  ))
  
  print(paste("Common genes across all matrices:", length(common_genes)))
  
  # Check if reasonable number of common genes
  if (length(common_genes) < 100) {
    print("WARNING: Few common genes across matrices, using pairwise intersection")
    # Try pairwise intersection instead
    common_bulk_ct <- intersect(rownames(bulk), rownames(cellTypeExpr))
    common_sc_ct <- intersect(rownames(singleCellExpr), rownames(cellTypeExpr))
    
    # Use whichever has more genes
    if (length(common_bulk_ct) > length(common_sc_ct)) {
      common_genes <- common_bulk_ct
    } else {
      common_genes <- common_sc_ct
    }
    print(paste("Using", length(common_genes), "common genes after pairwise intersection"))
  }
  
  # Sort for consistent ordering
  common_genes <- sort(common_genes)
  
  # Subset all matrices
  bulk <- bulk[common_genes, , drop=FALSE]
  singleCellExpr <- singleCellExpr[common_genes, , drop=FALSE]
  cellTypeExpr <- cellTypeExpr[common_genes, , drop=FALSE]
  
  # Force identical dimensions for feat lookup
  dimnames(bulk) <- list(common_genes, colnames(bulk))
  dimnames(cellTypeExpr) <- list(common_genes, colnames(cellTypeExpr))
  
  print("Matrix dimensions aligned for PREDE:")
  print(paste("bulk:", paste(dim(bulk), collapse="x")))
  print(paste("cellTypeExpr:", paste(dim(cellTypeExpr), collapse="x")))
}

# RNA-Sieve - Fix array dimension mismatch
if (method == "RNA-Sieve") {
  print("Enhanced RNA-Sieve fix")
  
  # The error suggests a major dimension mismatch between arrays
  # Create pre-aggregated data by cell type to avoid this
  
  # Filter genes
  common_genes <- intersect(rownames(bulk), rownames(singleCellExpr))
  print(paste("Common genes:", length(common_genes)))
  
  # Subset both matrices to common genes and sort them
  common_genes <- sort(common_genes)
  bulk <- bulk[common_genes, , drop=FALSE]
  singleCellExpr <- singleCellExpr[common_genes, , drop=FALSE]
  
  # Create pre-aggregated single-cell matrix by cell type
  # RNA-Sieve may be trying to auto-aggregate which causes dimension mismatch
  print("Creating pre-aggregated single-cell data for RNA-Sieve")
  
  # Create a specialized class that includes labeling info
  RNASieve_data <- list(
    bulk = bulk,
    sc = singleCellExpr,
    labels = singleCellLabels,
    common_genes = common_genes,
    cell_types = unique(singleCellLabels)
  )
  
  # Save this for RNA-Sieve to use
  assign("RNASieve_data", RNASieve_data, envir = .GlobalEnv)
  
  # Create and store full expression matrix with labels
  labeled_sc_expr <- singleCellExpr
  colnames(labeled_sc_expr) <- singleCellLabels
  assign("labeled_sc_expr", labeled_sc_expr, envir = .GlobalEnv)
  
  # Create a dense version of the matrix just in case
  bulk_dense <- as.matrix(bulk)
  sc_expr_dense <- as.matrix(singleCellExpr)
  assign("bulk_dense", bulk_dense, envir = .GlobalEnv)
  assign("sc_expr_dense", sc_expr_dense, envir = .GlobalEnv)
  
  print("RNA-Sieve inputs prepared")
}

# TOAST - Fix marker gene handling
if (method == "TOAST") {
  print("Enhanced TOAST fix")
  
  # TOAST expects marker genes in a specific format and has array dimension issues
  # Create more robust markers and cell-type expression
  
  if (is.null(markers) || length(markers) < length(all_cell_types)) {
    print("Creating robust markers for TOAST")
    robust_markers <- list()
    for (ct in all_cell_types) {
      # Find genes highly expressed in this cell type
      expr_in_ct <- avg_expr_by_celltype[, ct]
      expr_in_others <- rowMeans(avg_expr_by_celltype[, colnames(avg_expr_by_celltype) != ct, drop=FALSE])
      
      # Calculate ratio of expression
      ratio <- expr_in_ct / (expr_in_others + 1e-10)
      
      # Get top genes
      top_genes <- names(sort(ratio, decreasing=TRUE)[1:min(50, length(ratio))])
      robust_markers[[ct]] <- top_genes
    }
    markers <- robust_markers
  }
  
  # TOAST needs numMarker to be at least 2D
  # Let's ensure every cell type has at least 10 markers
  min_markers_per_type <- 10
  for (ct in names(markers)) {
    if (length(markers[[ct]]) < min_markers_per_type) {
      # Add more markers if needed
      missing_count <- min_markers_per_type - length(markers[[ct]])
      print(paste("Adding", missing_count, "more markers for", ct))
      
      # Find additional markers based on expression
      extra_markers <- setdiff(rownames(singleCellExpr), unlist(markers))
      if (length(extra_markers) >= missing_count) {
        markers[[ct]] <- c(markers[[ct]], sample(extra_markers, missing_count))
      } else {
        # Just duplicate existing markers in worst case
        markers[[ct]] <- c(markers[[ct]], rep(markers[[ct]], length.out=missing_count))
      }
    }
  }
  
  # Validate all markers exist in the data matrices
  final_markers <- list()
  for (ct in names(markers)) {
    valid_markers <- intersect(markers[[ct]], rownames(bulk))
    if (length(valid_markers) >= 5) {
      final_markers[[ct]] <- valid_markers
    } else {
      print(paste("Warning: Not enough valid markers for", ct))
    }
  }
  markers <- final_markers
  
  # Create a matrix version of the markers for TOAST
  markerList <- list()
  for (k in 1:length(markers)) {
    ct <- names(markers)[k]
    markerList[[k]] <- markers[[ct]]
  }
  names(markerList) <- names(markers)
  assign("markerList", markerList, envir = .GlobalEnv)
  
  print("TOAST markers prepared:")
  print(sapply(markers, length))
}

# DESeq2 and EMeth - Fix matrix dimensions
if (method %in% c("DESeq2", "EMeth")) {
  print(paste("Enhanced", method, "fix"))
  
  # Both methods require perfect alignment between bulk and cellTypeExpr
  if (is.null(cellTypeExpr)) {
    print("Creating cellTypeExpr from singleCellExpr")
    cellTypeExpr <- avg_expr_by_celltype
  }
  
  # Get common genes and strictly enforce same dimensions
  common_genes <- intersect(rownames(bulk), rownames(cellTypeExpr))
  print(paste("Common genes:", length(common_genes)))
  
  if (length(common_genes) > 0) {
    common_genes <- sort(common_genes)
    bulk <- bulk[common_genes, , drop=FALSE]
    cellTypeExpr <- cellTypeExpr[common_genes, , drop=FALSE]
    
    # Verify dimensions match
    if (nrow(bulk) != nrow(cellTypeExpr)) {
      stop(paste("Matrix dimensions still don't match after processing:", 
                 nrow(bulk), "vs", nrow(cellTypeExpr)))
    }
    
    print(paste("Matrices aligned with", nrow(bulk), "common genes"))
  } else {
    stop("No common genes found between bulk and cellTypeExpr")
  }
}

# Fixes for other methods mentioned as fixable
if (method %in% c("ARIC", "AutoGeneS")) {
  print(paste("Enhanced", method, "fix"))
  
  # These methods have issues with gene ID indexing
  # Clean and standardize all gene IDs across datasets
  
  # Convert all gene IDs to the same format regardless of source
  standard_genes <- intersect(rownames(bulk), rownames(singleCellExpr))
  print(paste("Common genes:", length(standard_genes)))
  
  if (length(standard_genes) > 100) {
    common_genes <- sort(standard_genes)
    bulk <- bulk[common_genes, , drop=FALSE]
    singleCellExpr <- singleCellExpr[common_genes, , drop=FALSE]
    
    # For AutoGeneS specifically - KeyError suggests it expects DataFrames with indices
    if (method == "AutoGeneS") {
      # Create named vectors that match gene IDs exactly
      gene_map <- data.frame(
        gene_id = common_genes,
        stringsAsFactors = FALSE
      )
      assign("gene_map", gene_map, envir = .GlobalEnv)
    }
    
    print("Common gene sets standardized for method")
  } else {
    stop("Insufficient common genes for method")
  }
}

if (method %in% c("BisqueMarker", "DSA", "MCPcounter")) {
  print(paste("Enhanced", method, "fix"))
  
  # These methods need well-defined marker genes
  # Ensure each cell type has sufficient markers
  
  if (is.null(markers) || any(sapply(markers, length) < 5)) {
    print("Enhancing marker genes for marker-based methods")
    
    # Already done in the general preprocessing section
    print("Using enhanced markers generated earlier")
    
    # Verify marker quality
    print("Marker counts per cell type:")
    print(sapply(markers, length))
  }
  
  # For MCPcounter, create a special version of the markers
  if (method == "MCPcounter") {
    # Create a flat gene list with cell type annotations
    mcp_markers <- data.frame(
      symbol = unlist(markers),
      cellType = rep(names(markers), sapply(markers, length)),
      stringsAsFactors = FALSE
    )
    assign("mcp_markers", mcp_markers, envir = .GlobalEnv)
    print("MCPcounter marker format prepared")
  }
}

if (method %in% c("Deblender", "DeCompress", "spatialDWLS", "DWLS")) {
  print(paste("Enhanced", method, "fix"))
  
  # These methods have issues with matrix dimensions or duplicate rownames
  
  # Fix duplicate rownames problem
  for (mat_name in c("bulk", "singleCellExpr", "cellTypeExpr", "signature")) {
    mat <- get(mat_name)
    if (!is.null(mat)) {
      if (any(duplicated(rownames(mat)))) {
        print(paste("Found", sum(duplicated(rownames(mat))), "duplicate rownames in", mat_name))
        rownames(mat) <- make.unique(rownames(mat))
        assign(mat_name, mat)
      }
    }
  }
  
  # For DWLS family, prepare coefficients structure for standard errors
  if (method %in% c("spatialDWLS", "DWLS")) {
    # Create a standard format with zero standard errors
    print("Preparing coefficient structure for DWLS/spatialDWLS")
    
    # Generate common genes
    common_genes <- intersect(rownames(bulk), rownames(singleCellExpr))
    
    # Prepare matrices
    bulk <- bulk[common_genes, , drop=FALSE]
    singleCellExpr <- singleCellExpr[common_genes, , drop=FALSE]
    
    # Force dense matrix format to avoid sparse matrix issues
    bulk <- as.matrix(bulk)
    singleCellExpr <- as.matrix(singleCellExpr)
    
    print("DWLS-specific preprocessing complete")
  }
}

if (method %in% c("debCAM", "EPIC")) {
  print(paste("Enhanced", method, "fix"))
  
  # debCAM - simplex initialization issues
  if (method == "debCAM") {
    # Need to ensure enough distinct data points for simplex construction
    print("Applying special fix for debCAM")
    
    # For debCAM, ensure data has sufficient dimensions and variance
    min_genes <- 100
    if (nrow(bulk) < min_genes) {
      print("WARNING: debCAM needs sufficient genes - adding pseudo-variation")
      
      # Create slightly noisy copies of existing genes to increase dimensions
      genes_needed <- min_genes - nrow(bulk)
      if (genes_needed > 0) {
        noise_factor <- 0.05
        
        # Create noisy duplicates
        extra_rows <- matrix(0, nrow=genes_needed, ncol=ncol(bulk))
        for (i in 1:genes_needed) {
          src_idx <- sample(1:nrow(bulk), 1)
          extra_rows[i,] <- bulk[src_idx,] * (1 + rnorm(ncol(bulk), 0, noise_factor))
        }
        
        # Add rownames
        rownames(extra_rows) <- paste0("gene_", 1:genes_needed)
        
        # Combine
        bulk <- rbind(bulk, extra_rows)
      }
    }
  }
  
  # EPIC - signature gene issues
  if (method == "EPIC") {
    # Need to ensure signature has sufficient genes
    print("Applying special fix for EPIC")
    
    # EPIC needs reference/signature gene profiles
    if (is.null(signature)) {
      print("Creating signature for EPIC from single-cell data")
      signature <- avg_expr_by_celltype
    }
    
    # Align gene sets
    common_genes <- intersect(rownames(bulk), rownames(signature))
    
    if (length(common_genes) < 100) {
      print("WARNING: Few common genes for EPIC signature - creating custom signature")
      
      # Create signature directly from single-cell data
      # with genes matching bulk exactly
      signature <- matrix(0, nrow=nrow(bulk), ncol=length(all_cell_types))
      rownames(signature) <- rownames(bulk)
      colnames(signature) <- all_cell_types
      
      # Use single-cell data to create signature where genes overlap
      overlap_genes <- intersect(rownames(bulk), rownames(singleCellExpr))
      
      # Create aggregate expression by cell type
      for (ct in all_cell_types) {
        cells_of_type <- singleCellLabels == ct
        if (sum(cells_of_type) > 0) {
          # For overlapping genes, use actual single-cell averages
          for (gene in overlap_genes) {
            gene_idx <- which(rownames(bulk) == gene)
            if (length(gene_idx) > 0) {
              cell_idxs <- which(cells_of_type)
              signature[gene_idx, ct] <- mean(singleCellExpr[gene, cell_idxs])
            }
          }
          
          # For non-overlapping genes, use random values following the distribution
          non_overlap <- setdiff(rownames(bulk), overlap_genes)
          for (gene in non_overlap) {
            gene_idx <- which(rownames(bulk) == gene)
            if (length(gene_idx) > 0) {
              # Sample from distribution of this cell type's expressions
              signature[gene_idx, ct] <- max(0, mean(bulk[gene_idx,]) * runif(1, 0.5, 1.5))
            }
          }
        }
      }
    } else {
      # Use common genes only
      bulk <- bulk[common_genes, , drop=FALSE]
      signature <- signature[common_genes, , drop=FALSE]
    }
    
    print(paste("Final signature dimensions for EPIC:", paste(dim(signature), collapse="x")))
  }
}

# Generate references
all_cell_types <- unique(singleCellLabels)
n_cell_types <- length(all_cell_types)

# Generate reference - happens after method-specific fixes
reference <- generateReference(singleCellExpr, singleCellLabels, type="signature")
signature <- reference$signature

#########################################################
# RUN DECONVOLUTION WITH ENHANCED ERROR HANDLING
#########################################################

print("Pre-processing complete, starting deconvolution...")

# Set environment variables for performance in various packages
Sys.setenv(NUMEXPR_MAX_THREADS="64")  # For DeconPeaker
Sys.setenv(OMP_NUM_THREADS="32")      # For general OpenMP settings
Sys.setenv(OPENBLAS_NUM_THREADS="32") # For OpenBLAS

# Wrap the deconvolution in error handling to get more specific error info
tryCatch({
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
  print("CHECK: Deconvolution completed successfully")
}, error = function(e) {
  print(paste("ERROR in deconvolution:", e$message))
  # Try to extract more detailed info
  print("Call stack:")
  print(traceback())
  # Try to see if there's a container error message
  tryCatch({
    last_logs <- system("tail -20 /tmp/singularity*.out", intern = TRUE)
    print("Last container logs:")
    print(last_logs)
  }, error = function(e2) {})
  
  # Re-throw the error to stop execution if deconvolution fails
  stop(paste("Deconvolution failed for method", method, ":", e$message))
})

# Post-process results to ensure valid data
# Normalize deconvolution results (absolute values and sum to 1 per sample)
if (!is.null(deconvolutionResult[[method]]$P)) {
  prop_matrix <- deconvolutionResult[[method]]$P
  
  # First check for and handle any NA values
  if (any(is.na(prop_matrix))) {
    print("Warning: NA values found in proportion matrix - replacing with zeros")
    prop_matrix[is.na(prop_matrix)] <- 0
  }
  
  # Process each sample (row) individually
  for (i in 1:nrow(prop_matrix)) {
    # Take absolute values to avoid negative proportions
    prop_matrix[i, ] <- abs(prop_matrix[i, ])
    
    # Normalize to sum to 1
    row_sum <- sum(prop_matrix[i, ])
    if (!is.na(row_sum) && row_sum > 0) {  # Avoid division by zero
      prop_matrix[i, ] <- prop_matrix[i, ] / row_sum
    } else {
      # If all values are zero or NA, set as all zeros (will show as error in plots)
      prop_matrix[i, ] <- rep(0, ncol(prop_matrix))
      print(paste("Warning: Sample", rownames(prop_matrix)[i], "had zero sum - keeping as zeros"))
    }
  }
  
  # Update the proportions matrix in the result object
  deconvolutionResult[[method]]$P <- prop_matrix
  print("Normalized cell type proportions (absolute values, sum to 1 per sample)")
} else {
  print("WARNING: No proportion matrix found in results!")
  stop(paste("Deconvolution method", method, "didn't produce proportion matrix (P)"))
}

#########################################################
# SAVE RESULTS
#########################################################

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