#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
  stop(paste("3 arguments must be supplied instead of", length(args)), call. = FALSE)
}

####### Parameter of script (ORDER IS IMPORTANT)
path_Rlibrary <- args[1] #IMPORTANT
input_data <- args[2]
output_data <- args[3]

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
# Save
output_file <- file.path(output_data, "bulk.rda")
save(bulk, file = output_file)
print(paste("Bulk matrix saved to:", output_file))