#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
  stop(paste("5 arguments must be supplied"), call. = FALSE)
}

# Parameters
path_Rlibrary <- args[1]  # R library path
seurat_file <- args[2]    # Path to Seurat RDS file
mapping_file <- args[3]   # Gene mapping file
output_dir <- args[4]     # Output directory
prefix <- args[5]         # Prefix for output files

# Libraries
.libPaths(path_Rlibrary, FALSE)
library(Seurat)

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
  cellTypeExpr[, i] <- rowMeans(ensembl_obj@assays$RNA@data[, cells_of_type]) # means vs sums plus normalize
}

# Save as RDA
rda_filename <- file.path(output_dir, paste0(prefix, "_celltype_expression.rda"))
save(cellTypeExpr, file=rda_filename)

# Save as CSV
csv_filename <- file.path(output_dir, paste0(prefix, "_celltype_expression.csv"))
cellTypeExpr_df <- as.data.frame(cellTypeExpr)
cellTypeExpr_df$EnsemblID <- rownames(cellTypeExpr_df)
write.csv(cellTypeExpr_df, file=csv_filename, row.names=FALSE)

print(paste("Cell type expression saved to:", rda_filename, "and", csv_filename))