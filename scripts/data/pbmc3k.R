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