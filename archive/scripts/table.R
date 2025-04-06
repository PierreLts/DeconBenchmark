#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript table.R <input_csv> <output_file> <benchmark_dir>", call.=FALSE)
}

# Parse command line arguments
input_csv <- args[1]  # Input CSV with model requirements
output_file <- args[2]  # Output file path
benchmark_dir <- args[3]  # Directory containing benchmark results

# Load required libraries
library(dplyr)
library(readr)
library(knitr)
library(stringr)
library(tools)

# Function to read benchmark results
read_benchmark_file <- function(file_path, suffix) {
  if (!file.exists(file_path)) {
    warning(paste("Benchmark file not found:", file_path))
    return(NULL)
  }
  
  cat(paste("Reading benchmark file:", file_path, "\n"))
  
  # Read benchmark data
  benchmark_data <- read_csv(file_path, show_col_types = FALSE)
  
  # Extract metrics of interest
  metrics <- c("JSD", "NRMSE", "PearsonCorr", "R2", "Runtime")
  
  # Find benchmark columns
  metric_cols <- c()
  model_cols <- c()
  
  for (metric in metrics) {
    # Look for columns that contain the metric name
    col_idx <- grep(paste0("^", metric, "$"), colnames(benchmark_data), ignore.case = TRUE)
    if (length(col_idx) > 0) {
      metric_cols <- c(metric_cols, col_idx[1])
    }
    
    # Look for "Models" columns that might have the metric name
    model_col_idx <- grep(paste0("^Models.*", metric, "$"), colnames(benchmark_data), ignore.case = TRUE)
    if (length(model_col_idx) > 0) {
      model_cols <- c(model_cols, model_col_idx[1])
    }
  }
  
  # Setup empty results dataframe
  result <- data.frame(Model = character(), stringsAsFactors = FALSE)
  
  # If less than 2 rows, there's not enough data
  if (nrow(benchmark_data) < 2) {
    warning("Benchmark file has insufficient data")
    return(result)
  }
  
  # Extract model information
  for (i in 2:nrow(benchmark_data)) {
    if (i > nrow(benchmark_data) || is.na(benchmark_data[i, 1])) next
    
    model_name <- as.character(benchmark_data[i, 1])
    if (model_name == "") next
    
    # Create new row for this model
    model_row <- data.frame(Model = model_name, stringsAsFactors = FALSE)
    
    # Add metrics columns with values
    for (j in seq_along(metric_cols)) {
      if (j > length(metric_cols)) next
      metric_name <- colnames(benchmark_data)[metric_cols[j]]
      col_name <- ifelse(suffix == "", metric_name, paste0(metric_name, suffix))
      model_row[[col_name]] <- as.character(benchmark_data[i, metric_cols[j]])
    }
    
    result <- rbind(result, model_row)
  }
  
  return(result)
}

# Read input models CSV file
cat("Reading input CSV file:", input_csv, "\n")
models_df <- read_csv(input_csv, show_col_types = FALSE)

# Read benchmark files
cat("Reading benchmark files...\n")
benchmark_files <- list(
  list(path = file.path(benchmark_dir, "TB_D100-bulk", "benchmarks", "TB_D100-bulk_benchmark_AB_select-AB.csv"), suffix = ""),
  list(path = file.path(benchmark_dir, "TB_D100-bulk_random", "benchmarks", "TB_D100-bulk_random_benchmark_AB_select-AB.csv"), suffix = "_Random"),
  list(path = file.path(benchmark_dir, "TB_D100-pseudobulk", "benchmarks", "TB_D100-pseudobulk_benchmark_AB_select-AB.csv"), suffix = "_Pseudo")
)

# Read each benchmark file and merge results
results_df <- models_df
for (file_info in benchmark_files) {
  cat(paste("Processing:", file_info$path, "\n"))
  benchmark_data <- read_benchmark_file(file_info$path, file_info$suffix)
  
  if (!is.null(benchmark_data) && nrow(benchmark_data) > 0) {
    # Merge with results (keep all rows even if not in benchmark)
    results_df <- left_join(results_df, benchmark_data, by = "Model")
  }
}

# Ensure numeric columns are formatted consistently
numeric_cols <- colnames(results_df)[sapply(results_df, is.numeric)]
for (col in numeric_cols) {
  results_df[[col]] <- as.character(results_df[[col]])
}

# Write the merged table to CSV
output_csv <- gsub("\\.html$", ".csv", output_file)
cat("Writing CSV output to:", output_csv, "\n")
write_csv(results_df, output_csv)

# Create a formatted HTML version if output is HTML
if (tolower(file_ext(output_file)) == "html") {
  cat("Generating HTML output...\n")
  
  # Create HTML table
  html_output <- kable(results_df, format = "html", 
                        caption = "Deconvolution Models Comparison")
  
  # Add some basic CSS styling
  html_style <- '
  <style>
    table {
      border-collapse: collapse;
      width: 100%;
      margin-bottom: 20px;
      font-family: Arial, sans-serif;
      font-size: 14px;
    }
    th, td {
      text-align: center;
      padding: 8px;
      border: 1px solid #ddd;
    }
    th {
      background-color: #4CAF50;
      color: white;
      font-weight: bold;
    }
    tr:nth-child(even) {
      background-color: #f2f2f2;
    }
    tr:hover {
      background-color: #ddd;
    }
    caption {
      font-size: 1.5em;
      margin-bottom: 10px;
      font-weight: bold;
      color: #333;
    }
    .requirement-1 {
      background-color: #d4edda;
      color: #155724;
      font-weight: bold;
    }
    .requirement-0 {
      background-color: #f8d7da;
      color: #721c24;
    }
    .model-name {
      text-align: left;
      font-weight: bold;
    }
    
    /* Performance metrics styling */
    .good-performance {
      background-color: #c3e6cb;
      font-weight: bold;
    }
    .medium-performance {
      background-color: #ffeeba;
    }
    .poor-performance {
      background-color: #f5c6cb;
    }
    
    /* Responsive design */
    @media screen and (max-width: 600px) {
      table {
        display: block;
        overflow-x: auto;
      }
    }
    
    /* Section headings */
    .section-heading {
      background-color: #343a40;
      color: white;
      font-weight: bold;
    }
  </style>
  '
  
  # Add color highlighting for requirements (1s and 0s)
  req_columns <- c("Ref-free", "Nb ct", "scExpr", "scLabels", "sc Subjects", "sc Subect Labels", 
                   "ctExpr", "Markers", "Signature", "Significant genes", "Methylation")
  
  for (col in req_columns) {
    html_output <- gsub(paste0(">1<"), paste0(" class=\"requirement-1\">1<"), html_output)
    html_output <- gsub(paste0(">0<"), paste0(" class=\"requirement-0\">0<"), html_output)
  }
  
  # Add explanatory text
  explanatory_text <- '
  <div style="margin-bottom: 20px;">
    <h2>Deconvolution Models Comparison Guide</h2>
    <p>This table helps you choose the right deconvolution model based on your available data and required performance.</p>
    
    <h3>Understanding the Requirements:</h3>
    <ul>
      <li><strong>Ref-free</strong>: 1 if the model is reference-free, 0 if it requires reference data</li>
      <li><strong>Nb ct</strong>: 1 if the model requires the number of cell types, 0 otherwise</li>
      <li><strong>scExpr</strong>: 1 if single-cell expression data is required, 0 otherwise</li>
      <li><strong>scLabels</strong>: 1 if single-cell labels are required, 0 otherwise</li>
      <li><strong>sc Subjects</strong>: 1 if single-cell subject information is required, 0 otherwise</li>
      <li><strong>ctExpr</strong>: 1 if cell type expression is required, 0 otherwise</li>
      <li><strong>Markers</strong>: 1 if marker genes are required, 0 otherwise</li>
      <li><strong>Signature</strong>: 1 if signature matrix is required, 0 otherwise</li>
    </ul>
    
    <h3>Performance Metrics:</h3>
    <ul>
      <li><strong>JSD</strong>: Jensen-Shannon Divergence (lower is better)</li>
      <li><strong>NRMSE</strong>: Normalized Root Mean Square Error (lower is better)</li>
      <li><strong>PearsonCorr</strong>: Pearson Correlation (higher is better)</li>
      <li><strong>R2</strong>: R-squared value (higher is better)</li>
      <li><strong>Runtime</strong>: Execution time in format HH:MM:SS</li>
    </ul>
    
    <p>Color coding: <span class="requirement-1">Required (1)</span> and <span class="requirement-0">Not Required (0)</span></p>
  </div>
  '
  
  html_full <- paste0("<!DOCTYPE html>\n<html>\n<head>\n", 
                      "<meta charset=\"UTF-8\">\n",
                      "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">\n",
                      "<title>Deconvolution Models Comparison</title>\n",
                      html_style, 
                      "</head>\n<body>\n",
                      explanatory_text,
                      html_output, 
                      "\n</body>\n</html>")
  
  # Write HTML output
  writeLines(html_full, output_file)
  cat("HTML output generated:", output_file, "\n")
}

cat("Table generation completed successfully!\n")