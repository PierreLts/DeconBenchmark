#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
  stop(paste("5 arguments must be supplied instead of", length(args)), call. = FALSE)
}

####### Parameter of script (ORDER IS IMPORTANT)
path_Rlibrary <- args[1] #IMPORTANT
input_data <- args[2]     # Path to PBMC Seurat object
output_dir <- args[3]     # Output directory
mapping_file <- args[4]   # Gene mapping file
prefix <- args[5]         # Prefix for output files

# Libraries
.libPaths(path_Rlibrary, FALSE) #IMPORTANT
library(Seurat)
library(Matrix)
print("CHECK: Libraries loaded")

# Load the PBMC dataset
print(paste("Loading Seurat object from:", input_data))
sc_data <- readRDS(input_data)
print("CHECK: Seurat object loaded")

# Print basic information about the object
print("=== Seurat Object Information ===")
print(sc_data)
print(paste("Number of cells:", ncol(sc_data)))
print(paste("Number of features:", nrow(sc_data)))
print(paste("Default assay:", DefaultAssay(sc_data)))

# Filter cells if needed (or use all cells)
filtered_sc <- sc_data

# Use BD_cell_type column for cell labels if available
if ("BD_cell_type" %in% colnames(filtered_sc@meta.data)) {
  print("Using BD_cell_type as cell type labels")
  cell_labels <- as.character(filtered_sc@meta.data$BD_cell_type)
} else if ("celltype" %in% colnames(filtered_sc@meta.data)) {
  print("Using celltype as cell type labels")
  cell_labels <- as.character(filtered_sc@meta.data$celltype)
} else {
  # Use active idents if no specific cell type column found
  print("No specific cell type column found, using active identities")
  cell_labels <- as.character(Idents(filtered_sc))
}

# Print cell type distribution
print("Cell type distribution:")
print(table(cell_labels, useNA = "ifany"))

# Create singleCellLabels
singleCellLabels <- cell_labels
names(singleCellLabels) <- colnames(filtered_sc)

# Create sample/subject information from orig.ident
print("Creating sample/subject information...")
if ("orig.ident" %in% colnames(filtered_sc@meta.data)) {
  print("Using 'orig.ident' as sample information")
  singleCellSubjects <- as.character(filtered_sc@meta.data$orig.ident)
  names(singleCellSubjects) <- colnames(filtered_sc)
} else {
  # Default sample name if orig.ident not available
  print("No sample information found. Using 'sample1' for all cells.")
  singleCellSubjects <- rep("sample1", length(singleCellLabels))
  names(singleCellSubjects) <- colnames(filtered_sc)
}

# Extract expression matrix
print("Extracting expression matrix...")
# Try to get counts data using Seurat's GetAssayData function
counts_data <- tryCatch({
  GetAssayData(filtered_sc, slot = "counts", assay = "RNA")
}, error = function(e) {
  # If that fails, try direct slot access
  tryCatch({
    filtered_sc@assays$RNA@counts
  }, error = function(e2) {
    # If that fails too, try accessing data in newer Seurat versions
    tryCatch({
      filtered_sc@assays$RNA@layers$counts
    }, error = function(e3) {
      stop("Could not extract count data using any method. Please check the Seurat object structure.")
    })
  })
})

# Convert to dense matrix if needed
if (inherits(counts_data, "dgCMatrix") || inherits(counts_data, "dgTMatrix")) {
  print("Converting sparse matrix to dense matrix...")
  singleCellExpr <- as.matrix(counts_data)
} else {
  singleCellExpr <- counts_data
}

print(paste("Expression matrix dimensions:", paste(dim(singleCellExpr), collapse=" x ")))

# Print first few gene names to understand their format
print("First 20 gene names in the dataset:")
head_genes <- head(rownames(singleCellExpr), 20)
print(head_genes)

# Load gene mapping file for ID conversion
if (file.exists(mapping_file)) {
  print(paste("Reading mapping file:", mapping_file))
  mapping_df <- read.table(mapping_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
  
  # Check mapping file structure
  print("Mapping file structure:")
  print(colnames(mapping_df))
  print(paste("Number of mappings:", nrow(mapping_df)))
  
  # Ensure proper column names for the first two columns
  if (ncol(mapping_df) >= 2) {
    colnames(mapping_df)[1:2] <- c("Ensembl_ID", "Gene_Name")
    
    # Remove duplicates from mapping
    mapping_df <- mapping_df[!duplicated(mapping_df$Gene_Name), ]
    print(paste("Unique mappings after removing duplicates:", nrow(mapping_df)))
    
    # Create a bidirectional mapping - both symbol->Ensembl and Ensembl->Ensembl
    gene_map <- setNames(mapping_df$Ensembl_ID, mapping_df$Gene_Name)
    
    # Add Ensembl->Ensembl identity mappings to ensure we keep existing Ensembl IDs
    ensembl_ids <- mapping_df$Ensembl_ID
    ensembl_self_map <- setNames(ensembl_ids, ensembl_ids)
    
    # Combine both mappings
    gene_map <- c(gene_map, ensembl_self_map)
    gene_map <- gene_map[!duplicated(names(gene_map))]
    
    print(paste("Final mapping entries:", length(gene_map)))
    
    # Apply gene mapping
    print("Applying sophisticated gene name mapping...")
    original_genes <- rownames(singleCellExpr)
    
    # Create a new mapping vector
    new_gene_names <- character(length(original_genes))
    
    # Process each gene
    ensembl_pattern <- "^ENSG\\d+"
    for (i in seq_along(original_genes)) {
      gene_name <- original_genes[i]
      
      # If it's already an Ensembl ID, keep it as is
      if (grepl(ensembl_pattern, gene_name)) {
        new_gene_names[i] <- gene_name
      } 
      # Otherwise, try to map it
      else if (gene_name %in% names(gene_map)) {
        new_gene_names[i] <- gene_map[gene_name]
      }
      # If mapping fails, keep the original
      else {
        new_gene_names[i] <- gene_name
      }
    }
    
    # Analyze mapping results
    already_ensembl <- sum(grepl(ensembl_pattern, original_genes))
    successfully_mapped <- sum(!grepl(ensembl_pattern, original_genes) & 
                              new_gene_names != original_genes)
    failed_to_map <- sum(!grepl(ensembl_pattern, original_genes) & 
                        new_gene_names == original_genes)
    
    print(paste("Mapping analysis:"))
    print(paste("- Already Ensembl IDs:", already_ensembl, 
                "(", round(already_ensembl/length(original_genes)*100, 1), "%)"))
    print(paste("- Successfully mapped to Ensembl:", successfully_mapped,
                "(", round(successfully_mapped/length(original_genes)*100, 1), "%)"))
    print(paste("- Failed to map:", failed_to_map,
                "(", round(failed_to_map/length(original_genes)*100, 1), "%)"))
    
    # Sample of mapping results
    print("Sample of mapping results:")
    sample_indices <- sample(seq_along(original_genes), min(10, length(original_genes)))
    for (i in sample_indices) {
      print(paste(original_genes[i], "->", new_gene_names[i]))
    }
    
    # Display patterns for unmapped genes to help diagnosis
    if (failed_to_map > 0) {
      unmapped_genes <- original_genes[!grepl(ensembl_pattern, original_genes) & 
                                      new_gene_names == original_genes]
      print(paste("First 20 unmapped genes:"))
      print(head(unmapped_genes, 20))
      
      # Analyze patterns in unmapped genes
      print("Analyzing patterns in unmapped genes:")
      pattern_counts <- list(
        "RP11-" = sum(grepl("^RP11-", unmapped_genes)),
        "LINC" = sum(grepl("^LINC", unmapped_genes)),
        "MIR" = sum(grepl("^MIR", unmapped_genes)),
        "AL" = sum(grepl("^AL\\d+", unmapped_genes)),
        "AC" = sum(grepl("^AC\\d+", unmapped_genes)),
        "AP" = sum(grepl("^AP\\d+", unmapped_genes)),
        "Numeric" = sum(grepl("^\\d+$", unmapped_genes)),
        "Cell barcodes" = sum(grepl("-\\d+$", unmapped_genes))
      )
      for (pattern_name in names(pattern_counts)) {
        print(paste("  -", pattern_name, ":", pattern_counts[[pattern_name]], "genes"))
      }
    }
    
    # Update gene names in the expression matrix
    rownames(singleCellExpr) <- new_gene_names
  } else {
    print("Warning: Mapping file doesn't have expected columns, skipping mapping")
  }
} else {
  print(paste("Warning: Mapping file not found:", mapping_file))
}

# Handle duplicate genes by summing their expression values
print("Checking for duplicate genes after mapping...")
gene_names <- rownames(singleCellExpr)
duplicate_count <- sum(duplicated(gene_names))

if (duplicate_count > 0) {
  print(paste("Found", duplicate_count, "duplicated gene names. Handling by summing expression values."))
  
  # Get unique gene names
  unique_genes <- unique(gene_names)
  
  # Create a new matrix for the results
  result <- matrix(0, nrow = length(unique_genes), ncol = ncol(singleCellExpr))
  rownames(result) <- unique_genes
  colnames(result) <- colnames(singleCellExpr)
  
  # Sum expression values for each unique gene
  for (i in seq_along(unique_genes)) {
    if (i %% 1000 == 0) print(paste("Processing gene", i, "of", length(unique_genes)))
    
    gene <- unique_genes[i]
    indices <- which(gene_names == gene)
    
    if (length(indices) > 0) {
      # Sum expression values across all instances of this gene
      result[i, ] <- colSums(singleCellExpr[indices, , drop = FALSE])
    }
  }
  
  # Replace the original matrix with the deduplicated one
  singleCellExpr <- result
  print(paste("After deduplication, expression matrix dimensions:", paste(dim(singleCellExpr), collapse=" x ")))
} else {
  print("No duplicate genes found after mapping.")
}

# Setup output directories
output_dir_base <- output_dir
output_dir_prefix <- file.path(output_dir_base, prefix)

# Create output directories if needed
if (!dir.exists(output_dir_base)) {
  dir.create(output_dir_base, recursive = TRUE)
}
if (!dir.exists(output_dir_prefix)) {
  dir.create(output_dir_prefix, recursive = TRUE)
}

# Save outputs - cell labels
print("Saving cell labels...")
save_path <- file.path(output_dir_prefix, paste0(prefix, "_singleCellLabels_AB.rda"))
save(singleCellLabels, file = save_path)

labels_df <- data.frame(
  CellBarcode = names(singleCellLabels),
  CellType = singleCellLabels,
  stringsAsFactors = FALSE
)
csv_path <- file.path(output_dir_prefix, paste0(prefix, "_singleCellLabels_AB.csv"))
write.csv(labels_df, file = csv_path, row.names = FALSE)

# Save outputs - expression matrix
print("Saving expression matrix...")
save_path <- file.path(output_dir_prefix, paste0(prefix, "_singleCellExpr_AB.rda"))
save(singleCellExpr, file = save_path)

# Save a small subset as CSV
max_rows <- min(50, nrow(singleCellExpr))
max_cols <- min(50, ncol(singleCellExpr))
small_expr <- singleCellExpr[1:max_rows, 1:max_cols]
csv_path <- file.path(output_dir_prefix, paste0(prefix, "_singleCellExpr_AB_50x50.csv"))
write.csv(small_expr, file = csv_path)

# Save subjects information
print("Saving cell subjects data...")
save_path <- file.path(output_dir_prefix, paste0(prefix, "_singleCellSubjects_AB.rda"))
save(singleCellSubjects, file = save_path)

subjects_df <- data.frame(
  CellBarcode = names(singleCellSubjects),
  Subject = singleCellSubjects,
  stringsAsFactors = FALSE
)
csv_path <- file.path(output_dir_prefix, paste0(prefix, "_singleCellSubjects_AB.csv"))
write.csv(subjects_df, file = csv_path, row.names = FALSE)

print("PBMC data processing completed successfully!")