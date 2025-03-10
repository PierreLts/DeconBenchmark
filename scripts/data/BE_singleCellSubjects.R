#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
  stop("3 arguments must be supplied: R_LIBRARY_PATH OUTPUT_DIR PREFIX", call.=FALSE)
}

# Parameters
path_Rlibrary <- args[1]  # R library path
output_dir <- args[2]     # Output directory
prefix <- args[3]         # Prefix for file names

# Set library path
.libPaths(path_Rlibrary, FALSE)

# Load DeconBenchmark package
library(DeconBenchmark)
print("CHECK: DeconBenchmark loaded")

# Load BloodExample data
data(BloodExample)
print("CHECK: BloodExample data loaded")
print(paste("Using prefix:", prefix))

# Extract single cell labels and their names
cell_labels <- BloodExample$singleCellLabels
num_cells <- length(cell_labels)

# Check if cell labels have names
print(paste("Cell labels have names:", !is.null(names(cell_labels))))
if (is.null(names(cell_labels))) {
  # If no names, create cell identifiers
  cell_ids <- paste0("cell_", 1:num_cells)
} else {
  cell_ids <- names(cell_labels)
}

# Create a simple subject assignment based on cell types
unique_cell_types <- unique(cell_labels)
num_subjects <- min(length(unique_cell_types), 5)  # Use up to 5 subjects

# Create mapping of cell types to subjects
cell_to_subject <- setNames(
  paste0("subject", 1 + (seq_along(unique_cell_types) - 1) %% num_subjects),
  unique_cell_types
)

# Assign subjects based on cell type
singleCellSubjects <- cell_to_subject[cell_labels]

# Save the result
save(singleCellSubjects, file = file.path(output_dir, paste0(prefix, "_singleCellSubjects.rda")))

# Save as CSV with proper cell identifiers
subjects_df <- data.frame(
  Cell = cell_ids,
  Subject = singleCellSubjects,
  stringsAsFactors = FALSE
)
write.csv(subjects_df, 
          file = file.path(output_dir, paste0(prefix, "_singleCellSubjects.csv")),
          row.names = FALSE)

print(paste("Single cell subjects saved to:", file.path(output_dir, paste0(prefix, "_singleCellSubjects.rda"))))
print(paste("Single cell subjects CSV saved to:", file.path(output_dir, paste0(prefix, "_singleCellSubjects.csv"))))
print(paste("Number of subjects created:", length(unique(singleCellSubjects))))
print(paste("Subject distribution:"))
print(table(singleCellSubjects))