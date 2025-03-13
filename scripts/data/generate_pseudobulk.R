## This script generates a pseudobulked raw counts matrix from scRNAseq data.

## Load libraries
library(dplyr)
library(Seurat)
# BiocManager::install("sva")
library(sva)
library(readxl)


## Integrated dataset
# Load data
seu_int <- readRDS("/Users/stang/Desktop/R_analysis/single_cell/seurat_integrated.rds")

# Filter your cells
# This has already been done in this object

# Extract raw counts
cell_barcodes <- colnames(seu_int@assays$RNA@counts)
metadata <- seu_int@meta.data

# Ensure metadata has a column for cell barcodes (if not already present)
metadata$Cell_Barcode <- rownames(metadata)
# Extract raw counts from the Seurat object
counts_matrix <- as.matrix(seu_int@assays$RNA@counts)

# Ensure metadata has grouping information
grouping_variable <- "Sample_Name"  # Replace with your grouping variable
metadata <- metadata %>% dplyr::mutate(Sample_Name = metadata[[grouping_variable]])
# Add Sample_Name to each column of counts matrix
metadata <- metadata[match(colnames(counts_matrix), metadata$Cell_Barcode), ]
metadata <- metadata %>%
  filter(Sample_Name != "Undetermined", Sample_Name != "Multiplet")

# Pseudobulk by grouping cells by Sample_Name
pseudobulk_counts <- t(as.matrix(counts_matrix)) %>%
  as.data.frame() %>%
  mutate(Sample_Name = metadata$Sample_Name) %>%
  group_by(Sample_Name) %>%
  summarise(across(everything(), sum)) %>%
  column_to_rownames("Sample_Name") %>%
  as.matrix()
pseudobulk_counts <- readRDS("/Users/stang/Desktop/R_analysis/single_cell/pseudobulk_counts_matrix_no_batch_correction.rds")
pseudobulk_counts <- t(pseudobulk_counts)
# Remove the last two columns
pseudobulk_counts <- pseudobulk_counts[, -((ncol(pseudobulk_counts) - 1):ncol(pseudobulk_counts))]

# Correct batch effects with sva

# Define batch variable (e.g., sequencing batch)
batch <- metadata$orig.ident[match(colnames(pseudobulk_counts), metadata$Sample_Name)]

# Apply ComBat to correct batch effects
pseudobulk_corrected <- ComBat_seq(pseudobulk_counts, batch=batch)

# Output: batch-corrected raw counts
head(pseudobulk_corrected)

# Save output
write.csv(pseudobulk_corrected, "/Users/stang/Desktop/R_analysis/single_cell/pseudobulk_counts_120k.csv", row.names = TRUE)






# ## Batch1/downsampled dataset
# # Load data
# seu_test <- readRDS("/Users/stang/Desktop/R_analysis/single_cell/batch1.rds") ## Replace with your downsampled object

# # Filter your cells
# seu_test <- Seurat::PercentageFeatureSet(seu_test, 
#                                          pattern = "^MT-", 
#                                          col.name = "percent.mito")
# seu_test <- Seurat::PercentageFeatureSet(seu_test, 
#                                          pattern = "^RP[SL]",
#                                          col.name = "percent.ribo")
# seu_test <- Seurat::PercentageFeatureSet(seu_test,
#                                          pattern = "^HB[^(P)]",
#                                          col.name = "percent.globin")
# seu_test <- subset(seu_test, subset = nFeature_RNA > 200 & 
#                      nFeature_RNA < 5000 &
#                      percent.mito < 25) # mito filter is much higher than standard 10x chromium (8%). This is standard for BD Rhapsody

# ## Extract raw counts
# cell_barcodes <- colnames(seu_test@assays$RNA@counts)
# metadata <- seu_test@meta.data

# # Ensure metadata has a column for cell barcodes (if not already present)
# metadata$Cell_Barcode <- rownames(metadata)
# # Extract raw counts from the Seurat object
# counts_matrix <- as.matrix(seu_test@assays$RNA@counts)

# # Ensure metadata has grouping information
# grouping_variable <- "Sample_Name"  # Replace with your grouping variable
# metadata <- metadata %>% dplyr::mutate(Sample_Name = metadata[[grouping_variable]])
# # Add Sample_Name to each column of counts matrix
# metadata <- metadata[match(colnames(counts_matrix), metadata$Cell_Barcode), ]

# # Pseudobulk by grouping cells by Sample_Name
# pseudobulk_counts <- t(as.matrix(counts_matrix)) %>%
#   as.data.frame() %>%
#   mutate(Sample_Name = metadata$Sample_Name) %>%
#   group_by(Sample_Name) %>%
#   summarise(across(everything(), sum)) %>%
#   column_to_rownames("Sample_Name") %>%
#   as.matrix()
# # pseudobulk_counts <- readRDS("/Users/stang/Desktop/R_analysis/single_cell/pseudobulk_counts_matrix_no_batch_correction.rds")
# pseudobulk_counts <- t(pseudobulk_counts)
# # Remove the last two columns
# pseudobulk_counts <- pseudobulk_counts[, -((ncol(pseudobulk_counts) - 1):ncol(pseudobulk_counts))]

# # Note: no need to correct for batch effects as there is none

# # Save output
# write.csv(pseudobulk_counts, "/Users/stang/Desktop/R_analysis/single_cell/pseudobulk_counts_30k.csv", row.names = TRUE)
