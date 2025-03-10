#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop(paste("4 arguments must be supplied"), call. = FALSE)
}

# Parameters
path_Rlibrary <- args[1]  # R library path
markers_rda <- args[2]    # Path to markers RDA file
output_dir <- args[3]     # Output directory
prefix <- args[4]         # Prefix for output files

# Libraries
.libPaths(path_Rlibrary, FALSE)

# Load markers
load(markers_rda)

# Create significant genes list from markers
significantGenes <- unique(unlist(markers))

# Save as RDA
rda_filename <- file.path(output_dir, paste0(prefix, "_significant_genes.rda"))
save(significantGenes, file=rda_filename)

# Save as CSV
csv_filename <- file.path(output_dir, paste0(prefix, "_significant_genes.csv"))
write.csv(data.frame(EnsemblID=significantGenes), file=csv_filename, row.names=FALSE)

print(paste("Significant genes saved to:", rda_filename, "and", csv_filename))