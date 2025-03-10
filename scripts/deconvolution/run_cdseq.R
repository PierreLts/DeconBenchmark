#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
  stop(paste("3 arguments must be supplied instead of", length(args)), call. = FALSE)
}

####### Parameter of script (ORDER IS IMPORTANT)
path_Rlibrary <- args[1] #IMPORTANT
input_dir <- args[2]     # Directory containing the generated data files
output_base_dir <- args[3]   # Base output directory for results

# Libraries
.libPaths(path_Rlibrary, FALSE) #IMPORTANT
library(CDSeq)
print("CHECK: Libraries loaded")

# Extract dataset prefix from the input directory name
dataset_prefix <- basename(input_dir)

# Create dataset-specific output directory
output_dir <- file.path(output_base_dir, dataset_prefix)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
print(paste("Results will be saved to:", output_dir))

# Load individual files from the directory
bulk_path <- file.path(input_dir, paste0(dataset_prefix, "_bulk.rda"))
singleCellLabels_path <- file.path(input_dir, paste0(dataset_prefix, "_singleCellLabels.rda"))
singleCellExpr_path <- file.path(input_dir, paste0(dataset_prefix, "_singleCellExpr.rda"))

if (!file.exists(bulk_path)) {
  stop(paste("Bulk data file not found at:", bulk_path))
}
if (!file.exists(singleCellLabels_path)) {
  stop(paste("Single cell labels file not found at:", singleCellLabels_path))
}
if (!file.exists(singleCellExpr_path)) {
  stop(paste("SingleCellExpr file not found at:", singleCellExpr_path))
}


# Load the data
load(bulk_path)
load(singleCellLabels_path)
load(singleCellExpr_path)

print("CHECK: Data loaded successfully")


# Get common genes between bulk and single-cell data
common_genes <- intersect(rownames(bulk), rownames(singleCellExpr))
message(paste("Common genes between datasets:", length(common_genes)))
# Subset both matrices to common genes
bulk <- bulk[common_genes, , drop=FALSE]
singleCellExpr <- singleCellExpr[common_genes, , drop=FALSE]


# Get number of cell types from the single-cell data
n_cell_types <- length(unique(singleCellLabels))
print(paste("Number of cell types detected:", n_cell_types))

# Run CDSeq
print("Starting CDSeq deconvolution...")
cdseq.result <- CDSeq(
  bulk_data = bulk,
  cell_type_number = n_cell_types, #seq(max(2, n_cell_types-3), min(n_cell_types+3, 15), 3),  # Try a range around the known number
  beta = 0.5,                 
  alpha = 5,                           
  mcmc_iterations = 700,             # (700-2000)
  dilution_factor = 5,               # (2-10)
  gene_subset_size = 300,            # (200-500)
  block_number = 5,                  # (>5)
  cpu_number = 16                    # CPUs
)

print("CHECK: CDSeq deconvolution completed, now formatting...")


# Structure the results to match plotting format
deconvolutionResult <- list()
deconvolutionResult$CDSeq <- list(
  P = t(cdseq.result$estProp),  # TRANSPOSED: Now samples in rows, cell types in columns
  S = cdseq.result$estGEP       # Genes in rows and cell types in columns
)

# Optionally, if you want to add annotation step



# Find common genes and subset both matrices
common_genes <- intersect(rownames(cdseq.result$estGEP), rownames(singleCellExpr))
print(paste("Number of common genes:", length(common_genes)))

# Subset both matrices to common genes
cdseq_gep_filtered <- cdseq.result$estGEP[common_genes, , drop=FALSE]
sc_gep_filtered <- singleCellExpr[common_genes, , drop=FALSE]

# Make sure the order is the same
sc_gep_filtered <- sc_gep_filtered[order(rownames(sc_gep_filtered)),]
cdseq_gep_filtered <- cdseq_gep_filtered[order(rownames(cdseq_gep_filtered)),]

# Check if the gene names are now matching
print("First 5 gene names from cdseq_gep_filtered:")
print(head(rownames(cdseq_gep_filtered), 5))
print("First 5 gene names from sc_gep_filtered:")
print(head(rownames(sc_gep_filtered), 5))


cdseq.result.celltypeassign <- cellTypeAssignSCRNA(
  cdseq_gep = cdseq_gep_filtered,
  cdseq_prop = cdseq.result$estProp,
  sc_gep = sc_gep_filtered,
  sc_annotation = singleCellLabels,
  sc_pt_size = 3,
  cdseq_pt_size = 6,
  seurat_nfeatures = 100,
  seurat_npcs = 50,
  seurat_dims=1:5,
  plot_umap = 0,
  plot_tsne = 0
)

# If successful, add annotated results to the output structure
if (!is.null(cdseq.result.celltypeassign) && !is.null(cdseq.result.celltypeassign$mapping)) {
  # Store the mapping
  deconvolutionResult$CDSeq$annotated_mapping <- cdseq.result.celltypeassign$mapping
  print("cell type assign succesful 1/2")
  # Apply the mapping to the column names of the transposed proportion matrix
  # This ensures the mapping aligns with the transposed format
  if (length(colnames(deconvolutionResult$CDSeq$P)) == length(cdseq.result.celltypeassign$mapping)) {
    # Replace numeric column names with actual cell type names
    new_colnames <- as.character(cdseq.result.celltypeassign$mapping)
    names(new_colnames) <- colnames(deconvolutionResult$CDSeq$P)
    colnames(deconvolutionResult$CDSeq$P) <- new_colnames
    print("cell type assign succesful 2/2")
  }
}


# Save results as RDA
results_filename <- file.path(output_dir, "results_CDSeq.rda")
save(deconvolutionResult, file=results_filename, compress=TRUE)

# Save proportions as CSV
proportions_df <- as.data.frame(deconvolutionResult$CDSeq$P)
proportions_df$Sample <- rownames(proportions_df)
csv_filename <- file.path(output_dir, "results_CDSeq_proportions.csv")
write.csv(proportions_df, file=csv_filename, row.names=FALSE)

print(paste("Results saved to:", results_filename))
print(paste("Proportions CSV saved to:", csv_filename))