#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
  stop(paste("5 arguments must be supplied"), call. = FALSE)
}

# Parameters
path_Rlibrary <- args[1]  # R library path
seurat_file <- args[2]    # Path to Seurat RDS file
output_dir <- args[3]     # Output directory
mapping_file <- args[4]   # Gene mapping file
prefix <- args[5]         # Prefix for output files

# Libraries
.libPaths(path_Rlibrary, FALSE)
library(Seurat)
library(dplyr)

# Load Seurat object
seurat_obj <- readRDS(seurat_file)

# Filter cells
valid_cells <- !seurat_obj@meta.data$Sample_Tag %in% c("Undetermined", "Multiplet")
filtered_obj <- subset(seurat_obj, cells = rownames(seurat_obj@meta.data)[valid_cells])

# Load gene mapping file
mapping_df <- read.table(mapping_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
colnames(mapping_df) <- c("Ensembl_ID", "Gene_Name")
mapping_df <- mapping_df[!duplicated(mapping_df$Gene_Name), ]

# Map gene names to Ensembl IDs
gene_names <- rownames(filtered_obj@assays$RNA@counts)
gene_map <- setNames(mapping_df$Ensembl_ID, mapping_df$Gene_Name)

# Find marker genes
Idents(filtered_obj) <- filtered_obj@meta.data$Cell_Type_Experimental
all_markers <- FindAllMarkers(filtered_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Create marker list
markers <- list()
markers_df_list <- list()

for (ct in unique(filtered_obj@meta.data$Cell_Type_Experimental)) {
  ct_markers <- all_markers %>% 
    filter(cluster == ct & p_val_adj < 0.05) %>%
    top_n(n = 50, wt = avg_log2FC) %>%
    pull(gene)
  
  # Skip if no markers found
  if (length(ct_markers) == 0) {
    print(paste("No markers found for cell type:", ct))
    next
  }
  
  # Map to Ensembl IDs and handle NAs
  ensembl_markers <- gene_map[ct_markers]
  valid_idx <- !is.na(ensembl_markers)
  
  if (sum(valid_idx) > 0) {
    # Only use markers with valid Ensembl IDs
    valid_symbols <- ct_markers[valid_idx]
    valid_ensembl <- ensembl_markers[valid_idx]
    
    markers[[ct]] <- valid_ensembl
    
    # Create data frame for this cell type
    temp_df <- data.frame(
      CellType = rep(ct, length(valid_symbols)),
      GeneSymbol = valid_symbols,
      EnsemblID = valid_ensembl,
      stringsAsFactors = FALSE
    )
    markers_df_list[[ct]] <- temp_df
  }
}

# Combine all marker data frames
if (length(markers_df_list) > 0) {
  markers_df <- do.call(rbind, markers_df_list)
}

# Convert to list format for deconvolution methods
markers_list <- list()
for (ct in unique(markers_df$CellType)) {
  ct_markers <- markers_df$EnsemblID[markers_df$CellType == ct]
  if (length(ct_markers) > 0) {
    markers_list[[ct]] <- ct_markers
  }
}

# Replace markers with properly structured list
markers <- markers_list

# Save markers as RDA
rda_filename <- file.path(output_dir, paste0(prefix, "_markers.rda"))
save(markers, file=rda_filename)

# Create a simplified CSV that matches the list structure
csv_list_format <- data.frame(
  CellType = rep(names(markers), sapply(markers, length)),
  MarkerID = unlist(markers)
)

# Save both CSV formats
csv_filename <- file.path(output_dir, paste0(prefix, "_markers.csv"))
write.csv(markers_df, file=csv_filename, row.names=FALSE)

# Save the list-like CSV format
csv_list_filename <- file.path(output_dir, paste0(prefix, "_markers_list.csv"))
write.csv(csv_list_format, file=csv_list_filename, row.names=FALSE)