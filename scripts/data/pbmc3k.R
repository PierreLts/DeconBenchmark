#!/usr/bin/Rscript
# Source data methodology: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

# Process command line arguments
args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3) {
  stop("At least 3 arguments must be supplied: R_LIBRARY_PATH INPUT_DATA_PATH OUTPUT_FILE_PATH", call.=FALSE)
}

# Set parameters from command line arguments
path_Rlibrary <- args[1]  # R library path
input_data_path <- args[2]  # Path to 10X data directory
output_file_path <- args[3]  # Path to save the output RDS file

# Set library path
.libPaths(path_Rlibrary, FALSE)

# Load required libraries
suppressPackageStartupMessages({
  library(biomaRt)
  library(dplyr)
  library(Seurat)
})

cat("Starting PBMC dataset processing\n")
cat("Input data path:", input_data_path, "\n")
cat("Output file path:", output_file_path, "\n")

# Load the PBMC dataset
cat("Reading 10X data...\n")
pbmc.data <- Read10X(data.dir = input_data_path)

# Initialize the Seurat object with the raw (non-normalized data)
cat("Creating Seurat object...\n")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
cat("Initial Seurat object dimensions:", dim(pbmc), "\n")

# Filtering steps
cat("Calculating mitochondrial gene percentage...\n")
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
cat("Filtering cells...\n")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
cat("After filtering, Seurat object dimensions:", dim(pbmc), "\n")

# Normalize data
cat("Normalizing data...\n")
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Find variable features
cat("Finding variable features...\n")
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Scale data
cat("Scaling data...\n")
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# Run PCA
cat("Running PCA...\n")
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Cluster cells
cat("Finding clusters...\n")
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Run UMAP
cat("Running UMAP...\n")
pbmc <- RunUMAP(pbmc, dims = 1:10)

# Annotate clusters according to tutorial
cat("Annotating cell types...\n")
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
pbmc$celltype <- Idents(pbmc)

# Convert metadata cell type labels to BD Rhapsody cell type labels
cat("Mapping to BD Rhapsody cell types...\n")
bd_mapping <- c(
  "T_CD4_naive",
  "Monocyte_classical",
  "T_CD4_memory",
  "B",
  "T_CD8",
  "Monocyte_nonclassical",
  "Natural_killer",
  "Dendritic",
  "Other"
)

names(bd_mapping) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, bd_mapping)
pbmc$BD_cell_type <- Idents(pbmc)

# Create output directory if it doesn't exist
output_dir <- dirname(output_file_path)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

# Save RDS
cat("Saving Seurat object to:", output_file_path, "\n")
saveRDS(pbmc, file = output_file_path)

cat("PBMC processing completed\n")

# Print final cell type counts
cat("Cell type counts:\n")
print(table(pbmc$BD_cell_type))


# Add this section at the end of pbmc3k.R after saving the RDS file

# Load the saved RDS file to inspect its structure
cat("\n====== PBMC RDS FILE STRUCTURE INSPECTION ======\n")
cat("Loading the saved RDS file to examine its structure...\n")
pbmc_loaded <- readRDS(output_file_path)

# Basic structure information
cat("\n1. OBJECT STRUCTURE:\n")
cat("Object class:", class(pbmc_loaded), "\n")
cat("Main slots in Seurat object:", paste(names(pbmc_loaded), collapse=", "), "\n")
cat("Dimensions (features x cells):", paste(dim(pbmc_loaded), collapse=" x "), "\n")

# Assay data tables
cat("\n2. ASSAY DATA TABLES:\n")
cat("Available assays:", paste(names(pbmc_loaded@assays), collapse=", "), "\n")
for (assay_name in names(pbmc_loaded@assays)) {
  assay_obj <- pbmc_loaded@assays[[assay_name]]
  cat("\n  ASSAY:", assay_name, "\n")
  cat("    Data tables in this assay:", paste(names(assay_obj), collapse=", "), "\n")
  
  # For each layer in the assay (counts, data, scale.data)
  for (layer_name in names(assay_obj)) {
    layer_data <- assay_obj[[layer_name]]
    cat("\n    TABLE:", layer_name, "\n")
    cat("      Class:", class(layer_data), "\n")
    cat("      Dimensions (genes x cells):", paste(dim(layer_data), collapse=" x "), "\n")
    
    # Print the first few row names (genes)
    cat("      First 5 genes:", paste(head(rownames(layer_data), 5), collapse=", "), "\n")
    
    # Print the first few column names (cells)
    cat("      First 5 cells:", paste(head(colnames(layer_data), 5), collapse=", "), "\n")
    
    # For sparse matrices, convert a small subset to dense for display
    if (inherits(layer_data, "dgCMatrix")) {
      cat("      Data preview (5x5) - converted from sparse matrix:\n")
      small_preview <- as.matrix(layer_data[1:min(5, nrow(layer_data)), 1:min(5, ncol(layer_data))])
    } else {
      cat("      Data preview (5x5):\n")
      small_preview <- layer_data[1:min(5, nrow(layer_data)), 1:min(5, ncol(layer_data))]
    }
    print(small_preview)
  }
}

# Metadata table
cat("\n3. METADATA TABLE:\n")
cat("Dimensions (cells x variables):", paste(dim(pbmc_loaded@meta.data), collapse=" x "), "\n")
cat("Column names:", paste(colnames(pbmc_loaded@meta.data), collapse=", "), "\n")
cat("First 5 rows of metadata:\n")
print(head(pbmc_loaded@meta.data, 5))

# Cell type information
cat("\n4. CELL TYPE TABLES:\n")
if ("celltype" %in% colnames(pbmc_loaded@meta.data)) {
  cat("Cell type table from 'celltype' column:\n")
  celltype_counts <- table(pbmc_loaded$celltype)
  print(as.data.frame(celltype_counts))
}

if ("BD_cell_type" %in% colnames(pbmc_loaded@meta.data)) {
  cat("\nBD Rhapsody cell type table from 'BD_cell_type' column:\n")
  bd_celltype_counts <- table(pbmc_loaded$BD_cell_type)
  print(as.data.frame(bd_celltype_counts))
}

# Commands list
cat("\n5. COMMANDS HISTORY:\n")
if (length(pbmc_loaded@commands) > 0) {
  cat("Commands executed on this object:\n")
  for (i in seq_along(pbmc_loaded@commands)) {
    cat("  ", i, ": ", names(pbmc_loaded@commands)[i], "\n", sep="")
  }
} else {
  cat("No commands history stored in this object\n")
}

cat("\n====== END OF RDS STRUCTURE INSPECTION ======\n")