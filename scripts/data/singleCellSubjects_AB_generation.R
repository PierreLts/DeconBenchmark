#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
  stop(paste("5 arguments must be supplied instead of", length(args)), call. = FALSE)
}

####### Parameter of script (ORDER IS IMPORTANT)
path_Rlibrary <- args[1] #IMPORTANT
input_data <- args[2]
output_dir <- args[3]
prefix <- args[4]
sample_filter <- args[5]  # "A", "B", or "A,B"

# Libraries
.libPaths(path_Rlibrary, FALSE) #IMPORTANT
library(Seurat)
print("CHECK: Libraries loaded")

# Load the single-cell dataset
sc_data <- readRDS(input_data)
print("CHECK: Seurat object loaded")

# Filter out "Undetermined" and "Multiplet" cells
cell_filter <- !sc_data@meta.data$Sample_Tag %in% c("Undetermined", "Multiplet")
filtered_sc <- subset(sc_data, cells = rownames(sc_data@meta.data)[cell_filter])

# Apply A/B filtering based on sample_filter parameter
filter_type <- tolower(sample_filter)
print(paste("Applying filter:", filter_type))

# Extract sample names from Sample_Name column
sample_names <- filtered_sc@meta.data$Sample_Name

# Create A/B filter based on the parameter
if (filter_type == "a") {
  ab_filter <- grepl("\\d+A$", sample_names)
  filter_label <- "A"
} else if (filter_type == "b") {
  ab_filter <- grepl("\\d+B$", sample_names)
  filter_label <- "B"
} else if (filter_type %in% c("a,b", "b,a", "ab", "both")) {
  # Keep both A and B samples
  ab_filter <- grepl("\\d+[AB]$", sample_names)
  filter_label <- "AB"
} else {
  stop(paste("Invalid filter type:", filter_type, "- Must be 'A', 'B', or 'A,B'"))
}

# Apply the filter
ab_filtered_sc <- subset(filtered_sc, cells = rownames(filtered_sc@meta.data)[ab_filter])

# Get count of cells for each group
a_count <- sum(grepl("\\d+A$", ab_filtered_sc@meta.data$Sample_Name))
b_count <- sum(grepl("\\d+B$", ab_filtered_sc@meta.data$Sample_Name))
total_count <- ncol(ab_filtered_sc)
print(paste("Selected cells count - A:", a_count, "B:", b_count, "Total:", total_count))

# Extract subject information
singleCellSubjects <- as.character(ab_filtered_sc@meta.data$Sample_Name)
names(singleCellSubjects) <- colnames(ab_filtered_sc)

# Count cells by subject
subject_counts <- table(singleCellSubjects)
print("Subject distribution:")
print(subject_counts)

# Save as RDA with filter label
rda_filename <- file.path(output_dir, paste0(prefix, "_singleCellSubjects_", filter_label, ".rda"))
save(singleCellSubjects, file = rda_filename)

# Save cell-subject mapping as CSV
cell_barcodes <- colnames(ab_filtered_sc)
subjects_df <- data.frame(
  Cell = cell_barcodes, 
  Subject = singleCellSubjects,
  stringsAsFactors = FALSE
)
csv_filename <- file.path(output_dir, paste0(prefix, "_singleCellSubjects_", filter_label, ".csv"))
write.csv(subjects_df, file = csv_filename, row.names = FALSE)

print(paste("Filtered cell subjects saved to:", rda_filename, "and", csv_filename))
print(paste("Total filtered cell count:", length(singleCellSubjects)))
print(paste("Unique subjects:", length(unique(singleCellSubjects))))