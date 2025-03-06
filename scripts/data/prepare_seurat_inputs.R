#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop(paste("4 arguments must be supplied"), call. = FALSE)
}

####### Parameters
path_Rlibrary <- args[1]
seurat_file <- args[2]
output_dir <- args[3]
mapping_file <- args[4]

# Libraries
.libPaths(path_Rlibrary, FALSE)
library(Seurat)
library(dplyr)

# Load Seurat object
seurat_obj <- readRDS(seurat_file)

# Filter cells
valid_cells <- !seurat_obj@meta.data$Sample_Tag %in% c("Undetermined", "Multiplet")
filtered_obj <- subset(seurat_obj, cells = rownames(seurat_obj@meta.data)[valid_cells])

### Mapping
# Load gene mapping file
mapping_df <- read.table(mapping_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
colnames(mapping_df) <- c("Ensembl_ID", "Gene_Name")
mapping_df <- mapping_df[!duplicated(mapping_df$Gene_Name), ]

# Map gene names to Ensembl IDs
gene_names <- rownames(filtered_obj@assays$RNA@counts)
gene_map <- setNames(mapping_df$Ensembl_ID, mapping_df$Gene_Name)

# Extract subject information
singleCellSubjects <- as.character(filtered_obj@meta.data$Sample_Name)
save(singleCellSubjects, file=file.path(output_dir, "single_cell_subjects.rda"))

# Save cell-subject mapping
cell_barcodes <- colnames(filtered_obj)
subjects_df <- data.frame(Cell=cell_barcodes, Subject=singleCellSubjects)
write.csv(subjects_df, file=file.path(output_dir, "single_cell_subjects.csv"), row.names=FALSE)
###

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
  # Save as CSV 
  write.csv(markers_df, file=file.path(output_dir, "markers.csv"), row.names=FALSE)
}

# Save markers as RDA
save(markers, file=file.path(output_dir, "markers.rda"))

# Create significant genes list from markers
significantGenes <- unique(unlist(markers))
save(significantGenes, file=file.path(output_dir, "significant_genes.rda"))
write.csv(data.frame(EnsemblID=significantGenes), 
          file=file.path(output_dir, "significant_genes.csv"), row.names=FALSE)

# Map all gene names to Ensembl IDs for cell type expression matrix
new_gene_names <- gene_map[gene_names]
new_gene_names[is.na(new_gene_names)] <- gene_names[is.na(new_gene_names)]

# Create Seurat object with Ensembl IDs
ensembl_obj <- filtered_obj
rownames(ensembl_obj@assays$RNA@counts) <- unname(new_gene_names)
rownames(ensembl_obj@assays$RNA@data) <- unname(new_gene_names)

# Generate cell type expression matrix
genes <- rownames(ensembl_obj@assays$RNA@data)
cell_types <- unique(ensembl_obj@meta.data$Cell_Type_Experimental)
cellTypeExpr <- matrix(0, nrow=length(genes), ncol=length(cell_types))
rownames(cellTypeExpr) <- genes
colnames(cellTypeExpr) <- cell_types

for (i in 1:length(cell_types)) {
  ct <- cell_types[i]
  cells_of_type <- which(ensembl_obj@meta.data$Cell_Type_Experimental == ct)
  cellTypeExpr[, i] <- rowMeans(ensembl_obj@assays$RNA@data[, cells_of_type])
}

save(cellTypeExpr, file=file.path(output_dir, "celltype_expression.rda"))

# Save CSV version with mapping back to gene symbols
cellTypeExpr_df <- as.data.frame(cellTypeExpr)
cellTypeExpr_df$EnsemblID <- rownames(cellTypeExpr_df)
write.csv(cellTypeExpr_df, file=file.path(output_dir, "celltype_expression.csv"), row.names=FALSE)