# Source data and scripts: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

library(biomaRt)
library(dplyr)
library(Seurat)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "/work/gr-fe/lorthiois/DeconBenchmark/scripts/data/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

# Filtering steps (I don't know if the reference needs to be filtered like normal scRNAseq - Pierre to confirm)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalise data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Find variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# Run PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Cluster cells
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Run UMAP
pbmc <- RunUMAP(pbmc, dims = 1:10)

# Annotate clusters according to tutorial
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
pbmc$celltype <- Idents(pbmc)

# Convert metadata cell type labels to BD Rhapsody cell type labels (for comparison of ground truth)
# Note that we do not have just a CD8 T Cell category in BD Rhapsody, just T_CD8_memory of T_CD8_naive.
# Thus, for the cells labeled as CD8 T in this reference, I will just change the name to T_CD8. 
# @Pierre you will just have to merge T_CD8_memory of T_CD8_naive when comparing against ground truth with this reference.
# Note that there is no category in this PBMC dataset for T_gamma_delta

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


# Save RDS
saveRDS(pbmc, file = "/Users/stang/Desktop/R_analysis/Deconvolution/filtered_gene_bc_matrices/hg19/pbmc_reference.rds")