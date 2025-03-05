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

# Load gene mapping file
mapping_df <- read.table(mapping_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
colnames(mapping_df) <- c("Ensembl_ID", "Gene_Name")
mapping_df <- mapping_df[!duplicated(mapping_df$Gene_Name), ]

# Map gene names to Ensembl IDs
gene_names <- rownames(filtered_obj@assays$RNA@counts)
gene_map <- setNames(mapping_df$Ensembl_ID, mapping_df$Gene_Name)
new_gene_names <- gene_map[gene_names]
new_gene_names[is.na(new_gene_names)] <- gene_names[is.na(new_gene_names)]

# Create a copy with Ensembl IDs
ensembl_obj <- filtered_obj
rownames(ensembl_obj@assays$RNA@counts) <- unname(new_gene_names)
rownames(ensembl_obj@assays$RNA@data) <- unname(new_gene_names)

# Extract subject information
# Fix for single cell subjects
singleCellSubjects <- as.character(ensembl_obj@meta.data$Sample_Name)
save(singleCellSubjects, file=file.path(output_dir, "single_cell_subjects.rda"))

# Get cell barcodes for CSV export
cell_barcodes <- colnames(ensembl_obj)
subjects_df <- data.frame(Cell=cell_barcodes, Subject=singleCellSubjects)
write.csv(subjects_df, file=file.path(output_dir, "single_cell_subjects.csv"), row.names=FALSE)


# Find marker genes
Idents(filtered_obj) <- filtered_obj@meta.data$Cell_Type_Experimental
all_markers <- FindAllMarkers(filtered_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Create marker list with Ensembl IDs
markers <- list()
markers_df <- data.frame()
for (ct in unique(filtered_obj@meta.data$Cell_Type_Experimental)) {
  ct_markers <- all_markers %>% 
    filter(cluster == ct & p_val_adj < 0.05) %>%
    top_n(n = 50, wt = avg_log2FC) %>%
    pull(gene)
  
  # Map to Ensembl IDs
  ensembl_markers <- gene_map[ct_markers]
  ensembl_markers[is.na(ensembl_markers)] <- ct_markers[is.na(ensembl_markers)]
  markers[[ct]] <- ensembl_markers
  
  # Add to dataframe for CSV export
  temp_df <- data.frame(
    CellType = ct,
    GeneSymbol = ct_markers,
    EnsemblID = ensembl_markers
  )
  markers_df <- rbind(markers_df, temp_df)
}
save(markers, file=file.path(output_dir, "markers.rda"))
write.csv(markers_df, file=file.path(output_dir, "markers.csv"), row.names=FALSE)

# Create significant genes with Ensembl IDs
significantGenes <- unique(unlist(markers))
save(significantGenes, file=file.path(output_dir, "significant_genes.rda"))
write.csv(data.frame(EnsemblID=significantGenes), 
          file=file.path(output_dir, "significant_genes.csv"), row.names=FALSE)

# Generate cell type expression matrix with Ensembl IDs
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

# Save as CSV (with gene names for readability)
cellTypeExpr_df <- as.data.frame(cellTypeExpr)
cellTypeExpr_df$EnsemblID <- rownames(cellTypeExpr_df)
# Map back to gene symbols for easy viewing
ensembl_to_symbol <- setNames(names(gene_map), gene_map)
cellTypeExpr_df$GeneSymbol <- ensembl_to_symbol[cellTypeExpr_df$EnsemblID]
cellTypeExpr_df$GeneSymbol[is.na(cellTypeExpr_df$GeneSymbol)] <- cellTypeExpr_df$EnsemblID[is.na(cellTypeExpr_df$GeneSymbol)]
write.csv(cellTypeExpr_df, file=file.path(output_dir, "celltype_expression.csv"), row.names=FALSE)