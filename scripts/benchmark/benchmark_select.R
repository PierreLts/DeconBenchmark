#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 6) {
  stop("6 arguments must be supplied: R_LIBRARY_PATH DATASET_PREFIX METRICS_CSV_PATH OUTPUT_DIR SAMPLE_FILTER SELECTION", call.=FALSE)
}

# Parameters
path_Rlibrary <- args[1]  # R library path
dataset_prefix <- args[2]  # Dataset prefix (e.g., TB1)
metrics_csv_path <- args[3]  # Path to the detailed metrics CSV
output_dir <- args[4]  # Output directory
sample_filter <- args[5]  # Sample filter used when generating data: A, B, or AB
selection <- args[6]  # Selection for analysis: A, B, or AB

# Set library path
.libPaths(path_Rlibrary, FALSE)

# Load required libraries
library(dplyr)
library(tidyr)
library(readr)
library(stringr)

# Function to extract runtime from log files
extract_runtime <- function(job_id) {
  log_file <- paste0("/scratch/lorthiois/logs/", job_id, ".e")
  if (!file.exists(log_file)) {
    return(NA)  # Return NA if log file doesn't exist
  }
  
  # Read the log file and look for runtime information
  log_content <- readLines(log_file)
  runtime_line <- grep("Runtime: .* seconds", log_content, value = TRUE)
  
  if (length(runtime_line) > 0) {
    # Extract the seconds value
    runtime_seconds <- as.numeric(gsub(".*Runtime: (\\d+) seconds.*", "\\1", runtime_line[1]))
    return(runtime_seconds)
  }
  
  return(NA)  # Return NA if runtime information not found
}

# Function to find and read model-job mapping file, selecting only the most recent job for each model
read_job_mapping <- function(dataset_prefix, sample_filter) {
  # Try different possible mapping file locations
  mapping_file_paths <- c(
    paste0("/work/gr-fe/lorthiois/DeconBenchmark/logs/", dataset_prefix, "_", sample_filter, "/model_job_mapping.txt"),
    paste0("/work/gr-fe/lorthiois/DeconBenchmark/logs/", dataset_prefix, "/model_job_mapping.txt")
  )
  
  mapping_file <- NULL
  for (path in mapping_file_paths) {
    if (file.exists(path)) {
      mapping_file <- path
      break
    }
  }
  
  if (is.null(mapping_file)) {
    warning("No mapping file found")
    return(data.frame(Method = character(), JobID = character(), stringsAsFactors = FALSE))
  }
  
  # Read the file content
  mapping_lines <- readLines(mapping_file)
  
  # Initialize a hash map to store the latest job ID for each method
  latest_jobs <- list()
  
  # Process mapping entries in order (so later entries override earlier ones)
  for (line in mapping_lines) {
    if (grepl(":", line) && !startsWith(line, "#")) {
      parts <- strsplit(line, ":")[[1]]
      if (length(parts) == 2) {
        method <- trimws(parts[1])
        job_id <- trimws(parts[2])
        
        # Store this job ID for the method, replacing any previous entry
        latest_jobs[[method]] <- job_id
      }
    }
  }
  
  # Convert the hash map to a data frame
  mapping_data <- data.frame(
    Method = names(latest_jobs),
    JobID = unlist(latest_jobs),
    stringsAsFactors = FALSE
  )
  
  print(paste("Found", nrow(mapping_data), "unique methods in job mapping"))
  
  return(mapping_data)
}

# Format runtime in a readable format (HH:MM:SS)
format_runtime <- function(seconds) {
  if (is.na(seconds)) return("N/A")
  
  # Ensure seconds is a numeric value
  seconds <- as.numeric(seconds)
  
  # Calculate hours, minutes, and seconds components
  hours <- floor(seconds / 3600)
  minutes <- floor((seconds %% 3600) / 60)
  remaining_seconds <- floor(seconds %% 60)  # Using floor to ensure integer values
  
  # Always include hours with zero-padding for all components
  return(sprintf("%02d:%02d:%02d", hours, minutes, remaining_seconds))
}

# Validate selection parameter
if (!selection %in% c("A", "B", "AB")) {
  stop("Selection parameter must be one of: 'A', 'B', or 'AB'")
}

# Read detailed metrics CSV
metrics_data <- read_csv(metrics_csv_path)
print(paste("Loaded metrics data with", nrow(metrics_data), "rows"))

# Get runtime data
print("Reading job mapping and extracting runtimes...")
mapping_data <- read_job_mapping(dataset_prefix, sample_filter)
print(paste("Found", nrow(mapping_data), "methods in job mapping"))

# Extract runtimes for each method
runtime_data <- data.frame(
  Method = character(),
  Runtime = numeric(),
  stringsAsFactors = FALSE
)

if (nrow(mapping_data) > 0) {
  for (i in 1:nrow(mapping_data)) {
    method <- mapping_data$Method[i]
    job_id <- mapping_data$JobID[i]
    runtime <- extract_runtime(job_id)
    
    runtime_data <- rbind(runtime_data, data.frame(
      Method = method,
      Runtime = runtime,
      stringsAsFactors = FALSE
    ))
  }
}

print(paste("Extracted runtime for", sum(!is.na(runtime_data$Runtime)), "methods"))

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

# Filter samples based on selection parameter
filtered_data <- metrics_data
if (selection == "A") {
  filtered_data <- metrics_data %>% filter(SampleGroup == "A")
  print(paste("Selected", nrow(filtered_data), "samples from group A"))
} else if (selection == "B") {
  filtered_data <- metrics_data %>% filter(SampleGroup == "B")
  print(paste("Selected", nrow(filtered_data), "samples from group B"))
} else {
  # selection == "AB" - keep all samples with valid group
  filtered_data <- metrics_data %>% filter(SampleGroup %in% c("A", "B"))
  print(paste("Selected", nrow(filtered_data), "samples from groups A and B"))
}

# Select relevant metrics
filtered_data <- filtered_data %>%
  select(Method, Sample, SampleGroup, PearsonCorr, NRMSE, R2, JSD)

# Calculate means and standard deviations by method
summary_stats <- filtered_data %>%
  group_by(Method) %>%
  summarize(
    PearsonCorr_mean = mean(PearsonCorr, na.rm = TRUE),
    PearsonCorr_sd = sd(PearsonCorr, na.rm = TRUE),
    NRMSE_mean = mean(NRMSE, na.rm = TRUE),
    NRMSE_sd = sd(NRMSE, na.rm = TRUE),
    R2_mean = mean(R2, na.rm = TRUE),
    R2_sd = sd(R2, na.rm = TRUE),
    JSD_mean = mean(JSD, na.rm = TRUE),
    JSD_sd = sd(JSD, na.rm = TRUE),
    Sample_count = n(),
    .groups = "drop"
  )

# *** MODIFIED SECTION: Sort all methods by JSD (lowest to highest) ***
# This is the key change - we sort all models by JSD once
sorted_by_jsd <- summary_stats %>% arrange(JSD_mean)
model_order <- sorted_by_jsd$Method

# Define metrics and create the result table structure
metrics <- c("PearsonCorr", "NRMSE", "R2", "JSD")
result_cols <- c("Model")  # Start with a single Model column

# Create column names for metric values and SDs
prefix <- paste0(selection, "_")
for (metric in metrics) {
  result_cols <- c(result_cols, 
                  paste0(prefix, metric),
                  paste0(prefix, "SD_", metric))
}

# Add runtime column
result_cols <- c(result_cols, paste0(prefix, "Runtime"))

# Create empty result dataframe with properly named columns
result_df <- data.frame(matrix(ncol = length(result_cols), nrow = 0))
colnames(result_df) <- result_cols

# Calculate overall metrics (mean values across all samples)
overall_metrics <- filtered_data %>%
  summarize(
    PearsonCorr_mean = mean(PearsonCorr, na.rm = TRUE),
    PearsonCorr_sd = sd(PearsonCorr, na.rm = TRUE),
    NRMSE_mean = mean(NRMSE, na.rm = TRUE),
    NRMSE_sd = sd(NRMSE, na.rm = TRUE),
    R2_mean = mean(R2, na.rm = TRUE),
    R2_sd = sd(R2, na.rm = TRUE),
    JSD_mean = mean(JSD, na.rm = TRUE),
    JSD_sd = sd(JSD, na.rm = TRUE)
  )

# Create summary row
summary_row <- data.frame(matrix(NA, nrow = 1, ncol = length(result_cols)))
colnames(summary_row) <- result_cols

# Set the Model column to "Overall"
summary_row[1, "Model"] <- paste0("Overall (", selection, ")")

# Fill in overall metrics for each measure type
for (metric in metrics) {
  value_col <- paste0(prefix, metric)
  sd_col_name <- paste0(prefix, "SD_", metric)
  
  # Add overall means and standard deviations
  summary_row[1, value_col] <- overall_metrics[[paste0(metric, "_mean")]]
  summary_row[1, sd_col_name] <- overall_metrics[[paste0(metric, "_sd")]]
}

# Add initial runtime info column (will update the runtime value later)
runtime_col <- paste0(prefix, "Runtime")
summary_row[1, runtime_col] <- "N/A" # Temporary, will be updated with average

# *** MODIFIED SECTION: Use JSD-sorted order for all metrics ***
# Initialize result data with the right number of rows
max_models <- nrow(sorted_by_jsd)
result_data <- as.data.frame(matrix(NA, nrow = max_models, ncol = length(result_cols)))
colnames(result_data) <- result_cols

# First, fill in the Model column with method names in JSD-sorted order
for (i in 1:nrow(sorted_by_jsd)) {
  if (i <= nrow(result_data)) {
    # Get the method at this position in the sorted list
    method <- sorted_by_jsd$Method[i]
    
    # Set the model name in the Model column
    result_data[i, "Model"] <- method
  }
}

# Process all metrics using the same JSD-based sorting
for (metric in metrics) {
  mean_col <- paste0(metric, "_mean")
  sd_col <- paste0(metric, "_sd")
  
  # Column names for this metric
  value_col <- paste0(prefix, metric)
  sd_col_name <- paste0(prefix, "SD_", metric)
  
  # Fill in data in JSD-sorted order
  for (i in 1:nrow(sorted_by_jsd)) {
    if (i <= nrow(result_data)) {
      # For all metrics, use the same row index to maintain consistency
      result_data[i, value_col] <- sorted_by_jsd[[mean_col]][i]
      result_data[i, sd_col_name] <- sorted_by_jsd[[sd_col]][i]
    }
  }
}

# *** MODIFIED SECTION: Use JSD-sorted order for runtime as well ***
# Add Runtime processing - using the same JSD-sorted order
# Join runtime data with the method list
if (nrow(runtime_data) > 0) {
  # Create a mapping of runtime data by method
  runtime_lookup <- setNames(runtime_data$Runtime, runtime_data$Method)
  
  # Fill runtime data in the same JSD-sorted order
  runtime_col <- paste0(prefix, "Runtime")
  
  for (i in 1:nrow(sorted_by_jsd)) {
    if (i <= nrow(result_data)) {
      method <- sorted_by_jsd$Method[i]
      
      # Look up runtime for this method
      runtime_value <- runtime_lookup[method]
      if (!is.na(runtime_value)) {
        result_data[i, runtime_col] <- format_runtime(runtime_value)
      } else {
        result_data[i, runtime_col] <- "N/A"
      }
    }
  }
  
  # Calculate average runtime for the Overall row
  avg_runtime <- mean(runtime_data$Runtime, na.rm = TRUE)
  summary_row[1, runtime_col] <- format_runtime(avg_runtime)
}

# Combine summary and data
final_result <- rbind(summary_row, result_data)

# Clean up column names for output
clean_colnames <- colnames(final_result)
clean_colnames <- gsub(paste0("^", prefix, "SD_"), "SD ", clean_colnames)
clean_colnames <- gsub(paste0("^", prefix, "Runtime"), "Runtime (H:M:S)", clean_colnames)
clean_colnames <- gsub(paste0("^", prefix), "", clean_colnames)
colnames(final_result) <- clean_colnames

# Write to CSV with selection in filename
output_file <- file.path(output_dir, paste0(dataset_prefix, "_benchmark_", sample_filter, "_select-", selection, ".csv"))
write.csv(final_result, output_file, row.names = FALSE)

print(paste("Benchmark results saved to:", output_file))