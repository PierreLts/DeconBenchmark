#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop(paste("4 arguments must be supplied instead of", length(args)), call. = FALSE)
}

####### Parameter of script (ORDER IS IMPORTANT)
path_Rlibrary <- args[1] #IMPORTANT
input_data <- args[2]
output_dir <- args[3]
prefix <- args[4]  # New: prefix for output files

# Libraries
.libPaths(path_Rlibrary, FALSE) #IMPORTANT
library(Matrix)

feature_counts <- read.csv(input_data, row.names = 1, check.names = FALSE)
# Gene version removal
rownames(feature_counts) <- sub("\\..*", "", rownames(feature_counts))
# Dense matrix
bulk <- as.matrix(feature_counts)

# Checks
print(bulk[1:5, 1:5])
print(class(bulk))

# Save as RDA
rda_filename <- file.path(output_dir, paste0(prefix, "_bulk.rda"))
save(bulk, file = rda_filename)

# Save as CSV
csv_filename <- file.path(output_dir, paste0(prefix, "_bulk.csv"))
write.csv(bulk, file = csv_filename)

print(paste("Bulk matrix saved to:", rda_filename, "and", csv_filename))