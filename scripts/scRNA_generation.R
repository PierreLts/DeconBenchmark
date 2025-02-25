#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop(paste("4 arguments must be supplied instead of", length(args)), call. = FALSE)
}


####### Parameter of script (ORDER IS IMPORTANT)
path_Rlibrary <- args[1] #IMPORTANT
input_data <- args[2]
output_data <- args[3]
mapping_file <- args[4]

###### MuSiC
# Libraries
.libPaths(path_Rlibrary, FALSE) #IMPORTANT
library(Seurat)
print("CHECK: Libraries")


# Load the single-cell dataset
sc_data <- readRDS(input_data)

# Filter out "Undetermined" and "Multiplet" cells
cell_filter <- !sc_data@meta.data$Sample_Tag %in% c("Undetermined", "Multiplet")
filtered_sc <- subset(sc_data, cells = rownames(sc_data@meta.data)[cell_filter])

# # Downsample to 5% of cells
# num_cells <- ncol(filtered_sc)
# selected_indices <- seq(1, num_cells, length.out = floor(0.05 * num_cells))
# selected_cells <- colnames(filtered_sc)[selected_indices]
# filtered_sc <- subset(filtered_sc, cells = selected_cells)

# Remove unused factor levels
filtered_sc@meta.data$Sample_Tag <- droplevels(filtered_sc@meta.data$Sample_Tag)

# Load gene mapping file
mapping_df <- read.table(mapping_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
colnames(mapping_df) <- c("Ensembl_ID", "Gene_Name")

# Remove duplicates from mapping
mapping_df <- mapping_df[!duplicated(mapping_df$Gene_Name), ]

# Map gene names to Ensembl IDs
gene_names <- rownames(filtered_sc@assays$RNA@counts)
gene_map <- setNames(mapping_df$Ensembl_ID, mapping_df$Gene_Name)
new_gene_names <- gene_map[gene_names]
new_gene_names[is.na(new_gene_names)] <- gene_names[is.na(new_gene_names)]

# Assign Ensembl IDs as rownames
rownames(filtered_sc@assays$RNA@counts) <- unname(new_gene_names)

# Convert sparse matrix to dense matrix
singleCellExpr <- as.matrix(filtered_sc@assays$RNA@counts)


# Checks
print("Preview of expression matrix (6x6):")
print(singleCellExpr[1:6, 1:6])
print(paste("Matrix dimensions:", paste(dim(singleCellExpr), collapse=" x ")))
print(paste("Matrix class:", class(singleCellExpr)))
# Save
output_file <- file.path(output_data, "singleCellExpr.rda")
save(singleCellExpr, file = output_file)
print(paste("Expression matrix saved to:", output_file))



