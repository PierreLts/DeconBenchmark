#!/usr/bin/Rscript

# This script generates radar plots from benchmark files
# It extracts metrics from benchmark CSV files and plots them in a radar chart

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript radar_compare.R <R_library_path> <output_file> <file1> <file2> ...")
}

# Set R library path
r_library_path <- args[1]
.libPaths(r_library_path, FALSE)

# Get output file path
output_file <- args[2]

# Get benchmark file paths
benchmark_files <- args[3:length(args)]

# Function to install packages if not available
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing package:", pkg, "\n")
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

# Install and load required packages
install_if_missing("fmsb")
install_if_missing("readr")
install_if_missing("dplyr")
install_if_missing("stringr")

# Function to extract dataset name and selection from file path
extract_name_and_selection <- function(file_path) {
  base_name <- basename(file_path)
  
  # Extract dataset name (everything before "_benchmark_")
  dataset_name <- gsub("_benchmark_.*$", "", base_name)
  
  # Extract selection (if present)
  if (grepl("select-", base_name)) {
    selection <- gsub(".*_select-([^\\.]+)\\.csv$", "\\1", base_name)
    return(paste0(dataset_name, "_select-", selection))
  }
  
  return(dataset_name)
}

# Read and process benchmark files
data_list <- list()
labels <- character()

for (file_path in benchmark_files) {
  # Extract dataset name and selection for labeling
  label <- extract_name_and_selection(file_path)
  labels <- c(labels, label)
  
  cat("Processing file:", file_path, "\n")
  
  # Read CSV file
  tryCatch({
    data <- read_csv(file_path)
    
    # Check file structure
    if (nrow(data) < 2) {
      warning("File has fewer than 2 rows: ", file_path)
      next
    }
    
    # Extract the overall values row (second row)
    metrics_row <- data[2, ]
    
    # Find columns for metrics
    # First try exact matches
    pearson_col <- which(colnames(data) %in% c("PearsonCorr", "Pearson"))[1]
    jsd_col <- which(colnames(data) == "JSD")[1]
    nrmse_col <- which(colnames(data) == "NRMSE")[1]
    r2_col <- which(colnames(data) %in% c("R2", "R^2"))[1]
    
    # If not found, try case-insensitive partial matches
    if (is.na(pearson_col)) pearson_col <- grep("pearson", colnames(data), ignore.case = TRUE)[1]
    if (is.na(jsd_col)) jsd_col <- grep("jsd", colnames(data), ignore.case = TRUE)[1]
    if (is.na(nrmse_col)) nrmse_col <- grep("nrmse", colnames(data), ignore.case = TRUE)[1]
    if (is.na(r2_col)) r2_col <- grep("r2|r\\^2", colnames(data), ignore.case = TRUE)[1]
    
    # Check if all metrics were found
    if (any(is.na(c(pearson_col, jsd_col, nrmse_col, r2_col)))) {
      warning("Could not find all metrics in file: ", file_path)
      print(colnames(data))
      next
    }
    
    # Extract metrics values
    pearson_val <- as.numeric(metrics_row[[pearson_col]])
    jsd_val <- as.numeric(metrics_row[[jsd_col]])
    nrmse_val <- as.numeric(metrics_row[[nrmse_col]])
    r2_val <- as.numeric(metrics_row[[r2_col]])
    
    cat("  Extracted metrics - Pearson:", pearson_val, "JSD:", jsd_val, "NRMSE:", nrmse_val, "R2:", r2_val, "\n")
    
    # Transform values to ensure best values are at the outer edge of radar plot
    # For JSD and NRMSE, lower is better, so we invert them
    metrics <- c(
      "PearsonCorr" = pearson_val,  # Range: -1 to 1, 1 is best (outer)
      "JSD" = 1 - jsd_val,          # Range: 0 to 1, 0 is best (outer) -> transform to 1 to 0
      "NRMSE" = -nrmse_val,         # Range: 0 to some+, 0 is best (outer) -> transform to negative
      "R2" = r2_val                 # Range: -5 to 1, 1 is best (outer)
    )
    
    data_list[[label]] <- metrics
    
  }, error = function(e) {
    warning("Error processing file ", file_path, ": ", e$message)
  })
}

if (length(data_list) == 0) {
  stop("No valid data found in any of the provided files")
}

# Create data frame for radar plot
radar_data <- as.data.frame(do.call(rbind, data_list))
rownames(radar_data) <- labels

# Define max and min values for each metric - fixed scales
max_values <- c(
  "PearsonCorr" = 1,     # Perfect positive correlation
  "JSD" = 1,             # Transformed, so 1 means JSD=0 (perfect match)
  "NRMSE" = 0,           # Transformed, so 0 means NRMSE=0 (perfect match)
  "R2" = 1               # Perfect R2 score
)

min_values <- c(
  "PearsonCorr" = -1,    # Perfect negative correlation
  "JSD" = 0,             # Transformed, so 0 means JSD=1 (worst)
  "NRMSE" = -5,          # Transformed, so -5 means NRMSE=5 (bad)
  "R2" = -5              # Very bad R2 score
)

# Add max and min rows required by fmsb
radar_data <- rbind(max_values, min_values, radar_data)

# Set up colors
colors_palette <- c("#1E90FF", "#FF69B4", "#32CD32", "#FFA500", "#9370DB", "#FF6347")
colors <- colors_palette[1:length(labels)]

# Create the radar chart
pdf(output_file, width = 10, height = 10)
par(mar = c(2, 2, 2, 2))

radarchart(
  radar_data,
  pcol = colors,
  pfcol = adjustcolor(colors, alpha.f = 0.3),
  plwd = 2,
  cglcol = "black",
  cglty = 1,
  axislabcol = "black",
  caxislabels = c("-5", "-4", "-1", "0.2", "1"),  # Fixed scale labels
  axistype = 1,
  calcex = 1.2,
  vlcex = 1.2,
  title = "Performance Metrics Comparison"
)

# Add a legend
legend(
  "bottomright",
  legend = labels,
  col = colors,
  lty = 1,
  lwd = 2,
  pch = 20,
  bty = "n"
)

dev.off()
cat("Radar plot saved to:", output_file, "\n")