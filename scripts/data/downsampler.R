#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
  stop("5 arguments must be supplied: R_LIBRARY_PATH INPUT_DIR OUTPUT_DIR PREFIX CELLS_PER_TYPE", call.=FALSE)
}

# Parameters
path_Rlibrary <- args[1]    # R library path
input_dir <- args[2]        # Source directory
output_dir <- args[3]       # Output directory
prefix <- args[4]           # Dataset prefix
cells_per_type <- as.numeric(args[5])  # Number of cells per cell type

# Set library path
.libPaths(path_Rlibrary, FALSE)

# First, find all the filter patterns that exist
message("Scanning for available filter patterns...")
files <- list.files(input_dir, pattern = paste0("^", prefix, "_singleCellExpr.*\\.rda$"))
filters <- c()

for (file in files) {
  # Extract the filter part from the filename
  filter_match <- sub(paste0("^", prefix, "_singleCellExpr"), "", file)
  filter_match <- sub("\\.rda$", "", filter_match)
  
  # If it starts with underscore, remove it
  if (startsWith(filter_match, "_")) {
    filter_match <- substring(filter_match, 2)
  }
  
  # Empty means no filter
  filters <- c(filters, filter_match)
}

if (length(filters) == 0) {
  stop("No singleCellExpr files found in the input directory")
}

message(paste("Found the following filter patterns:", paste(filters, collapse=", ")))

# Process each filter
for (filter in filters) {
  filter_suffix <- ifelse(filter == "", "", paste0("_", filter))
  
  # Construct file paths
  sc_expr_path <- file.path(input_dir, paste0(prefix, "_singleCellExpr", filter_suffix, ".rda"))
  sc_labels_path <- file.path(input_dir, paste0(prefix, "_singleCellLabels", filter_suffix, ".rda"))
  sc_subjects_path <- file.path(input_dir, paste0(prefix, "_singleCellSubjects", filter_suffix, ".rda"))
  
  # Basic file diagnostics
  message(paste("Checking files for filter:", filter))
  
  if (!file.exists(sc_expr_path)) {
    message(paste("Expression file not found:", sc_expr_path))
    next
  } else {
    # Basic file diagnostics
    file_info <- file.info(sc_expr_path)
    message(paste("Expression file size:", file_info$size, "bytes"))
    if (file_info$size == 0) {
      message("WARNING: Expression file is empty!")
      next
    }
  }
  
  if (!file.exists(sc_labels_path)) {
    message(paste("Labels file not found:", sc_labels_path))
    next
  } else {
    file_info <- file.info(sc_labels_path)
    message(paste("Labels file size:", file_info$size, "bytes"))
    if (file_info$size == 0) {
      message("WARNING: Labels file is empty!")
      next
    }
  }
  
  # Try to detect file problems
  tryCatch({
    # Try to load just a small part of the file
    con <- file(sc_expr_path, "rb")
    header <- readBin(con, "raw", 100)  # Read first 100 bytes
    close(con)
    message("Successfully read file header")
  }, error = function(e) {
    message(paste("Warning: File access test failed:", e$message))
  })
  
  # Load the data with error handling
  message(paste("Processing filter:", ifelse(filter == "", "none", filter)))
  
  # Load expression data
  expr_loaded <- FALSE
  tryCatch({
    load(sc_expr_path)
    if (exists("singleCellExpr")) {
      message(paste("Expression matrix dimensions:", paste(dim(singleCellExpr), collapse=" x ")))
      expr_loaded <- TRUE
    } else {
      message("Error: The variable 'singleCellExpr' was not found in the file")
    }
  }, error = function(e) {
    message(paste("Error loading expression file:", e$message))
    message(paste("Failed to load:", sc_expr_path))
  })
  
  if (!expr_loaded) {
    message(paste("Skipping filter", filter, "due to expression data loading failure"))
    next
  }
  
  # Load labels data
  labels_loaded <- FALSE
  tryCatch({
    load(sc_labels_path)
    if (exists("singleCellLabels")) {
      message(paste("Number of cell labels:", length(singleCellLabels)))
      labels_loaded <- TRUE
    } else {
      message("Error: The variable 'singleCellLabels' was not found in the file")
    }
  }, error = function(e) {
    message(paste("Error loading labels file:", e$message))
    message(paste("Failed to load:", sc_labels_path))
  })
  
  if (!labels_loaded) {
    message(paste("Skipping filter", filter, "due to labels data loading failure"))
    next
  }
  
  # Verify cell barcodes match between expression and labels
  expr_cells <- colnames(singleCellExpr)
  labels_cells <- names(singleCellLabels)
  
  if (!all(expr_cells %in% labels_cells) || !all(labels_cells %in% expr_cells)) {
    message("WARNING: Cell barcodes in expression and labels do not match perfectly")
    message(paste("Expression unique cells:", length(expr_cells)))
    message(paste("Labels unique cells:", length(labels_cells)))
    message(paste("Common cells:", length(intersect(expr_cells, labels_cells))))
  }
  
  # Load subjects if available
  has_subjects <- file.exists(sc_subjects_path)
  if (has_subjects) {
    subjects_loaded <- FALSE
    tryCatch({
      load(sc_subjects_path)
      if (exists("singleCellSubjects")) {
        message(paste("Number of cell subjects:", length(singleCellSubjects)))
        subjects_loaded <- TRUE
      } else {
        message("Error: The variable 'singleCellSubjects' was not found in the file")
        has_subjects <- FALSE
      }
    }, error = function(e) {
      message(paste("Error loading subjects file:", e$message))
      message("Will continue without subjects data")
      has_subjects <- FALSE
    })
  } else {
    message("No subjects file found, will not process subjects")
  }
  
  # Get original counts
  cell_types <- unique(singleCellLabels)
  message(paste("Found", length(cell_types), "cell types"))
  
  original_counts <- table(singleCellLabels)
  message("Original cell counts:")
  print(original_counts)
  
  # Identify cells to keep
  cells_to_keep <- c()
  
  # Sample cells from each cell type
  for (cell_type in cell_types) {
    cells_of_type <- names(singleCellLabels[singleCellLabels == cell_type])
    n_cells <- length(cells_of_type)
    
    # Take all cells if fewer than requested, otherwise sample
    if (n_cells <= cells_per_type) {
      message(paste("Cell type", cell_type, "has only", n_cells, "cells (keeping all)"))
      sampled_cells <- cells_of_type
    } else {
      message(paste("Downsampling cell type", cell_type, "from", n_cells, "to", cells_per_type, "cells"))
      set.seed(42 + which(cell_types == cell_type))  # Use consistent seed per cell type
      sampled_cells <- sample(cells_of_type, cells_per_type)
    }
    
    cells_to_keep <- c(cells_to_keep, sampled_cells)
  }
  
  # Subset data
  singleCellExpr_ds <- singleCellExpr[, cells_to_keep, drop=FALSE]
  singleCellLabels_ds <- singleCellLabels[cells_to_keep]
  
  if (has_subjects && subjects_loaded) {
    singleCellSubjects_ds <- singleCellSubjects[cells_to_keep]
  }
  
  # Verify downsampled data
  downsampled_counts <- table(singleCellLabels_ds)
  message("Downsampled cell counts:")
  print(downsampled_counts)
  
  # Create directory if it doesn't exist
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Save the results
  singleCellExpr <- singleCellExpr_ds
  singleCellLabels <- singleCellLabels_ds
  
  # Save expression matrix RDA
  tryCatch({
    save(singleCellExpr, file = file.path(output_dir, paste0(prefix, "_singleCellExpr", filter_suffix, ".rda")))
    
    # Also save as CSV (small subset due to potentially large size)
    max_rows <- min(1000, nrow(singleCellExpr))
    max_cols <- min(100, ncol(singleCellExpr))
    expr_subset <- singleCellExpr[1:max_rows, 1:max_cols, drop=FALSE]
    write.csv(expr_subset, file = file.path(output_dir, paste0(prefix, "_singleCellExpr", filter_suffix, "_subset.csv")))
    message("Saved expression data (RDA and CSV subset)")
  }, error = function(e) {
    message(paste("Error saving expression data:", e$message))
  })
  
  # Save labels RDA and CSV
  tryCatch({
    save(singleCellLabels, file = file.path(output_dir, paste0(prefix, "_singleCellLabels", filter_suffix, ".rda")))
    
    # Create a data frame for the labels and save as CSV
    labels_df <- data.frame(
      CellID = names(singleCellLabels),
      CellType = singleCellLabels,
      stringsAsFactors = FALSE
    )
    write.csv(labels_df, file = file.path(output_dir, paste0(prefix, "_singleCellLabels", filter_suffix, ".csv")), row.names = FALSE)
    message("Saved labels data (RDA and CSV)")
  }, error = function(e) {
    message(paste("Error saving labels data:", e$message))
  })
  
  # Save subjects if available
  if (has_subjects && subjects_loaded) {
    tryCatch({
      singleCellSubjects <- singleCellSubjects_ds
      save(singleCellSubjects, file = file.path(output_dir, paste0(prefix, "_singleCellSubjects", filter_suffix, ".rda")))
      
      # Create a data frame for the subjects and save as CSV
      subjects_df <- data.frame(
        CellID = names(singleCellSubjects),
        Subject = singleCellSubjects,
        stringsAsFactors = FALSE
      )
      write.csv(subjects_df, file = file.path(output_dir, paste0(prefix, "_singleCellSubjects", filter_suffix, ".csv")), row.names = FALSE)
      message("Saved subjects data (RDA and CSV)")
    }, error = function(e) {
      message(paste("Error saving subjects data:", e$message))
    })
  }
  
  message(paste("Saved downsampled files for filter:", ifelse(filter == "", "none", filter)))
}

message("Downsampling completed successfully")