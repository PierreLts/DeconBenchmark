#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop("4 arguments must be supplied: R_LIBRARY_PATH DATASET_PREFIX METRICS_CSV_PATH OUTPUT_DIR", call.=FALSE)
}

# Parameters
path_Rlibrary <- args[1]  # R library path
dataset_prefix <- args[2]  # Dataset prefix (e.g., TB1)
metrics_csv_path <- args[3]  # Path to the detailed metrics CSV
output_dir <- args[4]  # Output directory

# Set library path
.libPaths(path_Rlibrary, FALSE)

# Load required libraries
library(dplyr)
library(tidyr)
library(readr)
library(stringr)

# Read detailed metrics CSV
metrics_data <- read_csv(metrics_csv_path)
print(paste("Loaded metrics data with", nrow(metrics_data), "rows"))

# Function to determine if a sample is from group A or B
get_sample_group <- function(sample_name) {
  # Look for A or B at the end of a number sequence
  if (grepl("\\d+A", sample_name)) return("A")
  if (grepl("\\d+B", sample_name)) return("B")
  
  # If no clear pattern, just look for the last character
  last_char <- substr(sample_name, nchar(sample_name), nchar(sample_name))
  if (last_char %in% c("A", "B")) return(last_char)
  
  return(NA)  # Undetermined group
}

# Add sample group column
metrics_data$SampleGroup <- sapply(metrics_data$Sample, get_sample_group)
print(paste("Found", sum(metrics_data$SampleGroup == "A"), "samples in group A"))
print(paste("Found", sum(metrics_data$SampleGroup == "B"), "samples in group B"))

# Filter to keep only A and B groups and relevant metrics
filtered_data <- metrics_data %>%
  filter(SampleGroup %in% c("A", "B")) %>%
  select(Method, Sample, SampleGroup, PearsonCorr, SpearmanCorr, MAE, RMSE, R2, JSD)

# Calculate means and standard deviations by method and group
summary_stats <- filtered_data %>%
  group_by(Method, SampleGroup) %>%
  summarize(
    PearsonCorr_mean = mean(PearsonCorr, na.rm = TRUE),
    PearsonCorr_sd = sd(PearsonCorr, na.rm = TRUE),
    SpearmanCorr_mean = mean(SpearmanCorr, na.rm = TRUE),
    SpearmanCorr_sd = sd(SpearmanCorr, na.rm = TRUE),
    MAE_mean = mean(MAE, na.rm = TRUE),
    MAE_sd = sd(MAE, na.rm = TRUE),
    RMSE_mean = mean(RMSE, na.rm = TRUE),
    RMSE_sd = sd(RMSE, na.rm = TRUE),
    R2_mean = mean(R2, na.rm = TRUE),
    R2_sd = sd(R2, na.rm = TRUE),
    JSD_mean = mean(JSD, na.rm = TRUE),
    JSD_sd = sd(JSD, na.rm = TRUE),
    Sample_count = n(),
    .groups = "drop"
  )

# Define metrics and create the result table structure
metrics <- c("PearsonCorr", "SpearmanCorr", "MAE", "RMSE", "R2", "JSD")
result_cols <- c()

# Create column names with prefixes to avoid duplicates
for (group in c("A", "B")) {
  for (metric in metrics) {
    prefix <- paste0(group, "_")
    result_cols <- c(result_cols, 
                    paste0(prefix, "Models_", metric),
                    paste0(prefix, metric),
                    paste0(prefix, "SD_", metric))
  }
}

# Create empty result dataframe with properly named columns
result_df <- data.frame(matrix(ncol = length(result_cols), nrow = 0))
colnames(result_df) <- result_cols

# Add header row with "A samples" and "B Samples"
header_row <- data.frame(matrix(ncol = length(result_cols), nrow = 1))
colnames(header_row) <- result_cols

# Fill header row
a_cols <- grep("^A_", result_cols)
b_cols <- grep("^B_", result_cols)
header_row[1, a_cols] <- "A samples"
header_row[1, b_cols] <- "B Samples"

# Process each group
max_models <- max(summary_stats %>% group_by(SampleGroup) %>% summarize(n = n()) %>% pull(n))
result_data <- as.data.frame(matrix(NA, nrow = max_models, ncol = length(result_cols)))
colnames(result_data) <- result_cols

# Process each group and metric
for (group in c("A", "B")) {
  group_data <- summary_stats %>% filter(SampleGroup == group)
  
  for (metric in metrics) {
    mean_col <- paste0(metric, "_mean")
    sd_col <- paste0(metric, "_sd")
    
    # Sort in appropriate direction
    if (metric %in% c("PearsonCorr", "SpearmanCorr", "R2")) {
      sorted_data <- group_data %>% arrange(desc(!!sym(mean_col)))
    } else {
      sorted_data <- group_data %>% arrange(!!sym(mean_col))
    }
    
    # Fill in the appropriate columns
    model_col <- paste0(group, "_Models_", metric)
    value_col <- paste0(group, "_", metric)
    sd_col_name <- paste0(group, "_SD_", metric)
    
    # Fill in data
    for (i in 1:nrow(sorted_data)) {
      if (i <= nrow(result_data)) {
        result_data[i, model_col] <- sorted_data$Method[i]
        result_data[i, value_col] <- sorted_data[[mean_col]][i]
        result_data[i, sd_col_name] <- sorted_data[[sd_col]][i]
      }
    }
  }
}

# Combine header and data
final_result <- rbind(header_row, result_data)

# Clean up column names for output
clean_colnames <- gsub("^A_Models_", "Models ", colnames(final_result))
clean_colnames <- gsub("^A_SD_", "SD ", clean_colnames)
clean_colnames <- gsub("^B_Models_", "Models ", clean_colnames)
clean_colnames <- gsub("^B_SD_", "SD ", clean_colnames)
clean_colnames <- gsub("^A_", "", clean_colnames)
clean_colnames <- gsub("^B_", "", clean_colnames)
colnames(final_result) <- clean_colnames

# Write to CSV with sample filter in filename
output_file <- file.path(output_dir, paste0(dataset_prefix, "_benchmark_", sample_filter, ".csv"))
write.csv(final_result, output_file, row.names = FALSE)

print(paste("Benchmark results saved to:", output_file))