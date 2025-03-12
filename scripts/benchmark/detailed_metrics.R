#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
  stop(paste("5 arguments must be supplied instead of", length(args)), call. = FALSE)
}

####### Parameter of script
path_Rlibrary <- args[1] #IMPORTANT
dataset_prefix <- args[2] # Dataset prefix (e.g., TB1)
output_base_dir <- args[3] # Base output directory for benchmark results
include_overall_gt <- as.logical(args[4]) # Whether to include overall ground truth
sample_filter <- args[5] # Sample filter: A, B, or AB

# Libraries
.libPaths(path_Rlibrary, FALSE) #IMPORTANT
library(reshape2)
library(dplyr)
library(tidyr)
library(Metrics)  # For common error metrics
print("Check Libraries")

# Setup paths
results_dir <- file.path("/work/gr-fe/lorthiois/DeconBenchmark/deconv_results", dataset_prefix)
data_dir <- file.path("/work/gr-fe/lorthiois/DeconBenchmark/generated_data", dataset_prefix)
output_dir <- file.path(output_base_dir, dataset_prefix, "benchmarks")

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load per-sample ground truth without filter suffix
per_sample_gt_path <- file.path(data_dir, paste0(dataset_prefix, "_GT_proportions_per_sample.rda"))
if (!file.exists(per_sample_gt_path)) {
  stop(paste("Per-sample ground truth file not found:", per_sample_gt_path))
}
load(per_sample_gt_path)
per_sample_gt_proportions <- groundTruth$P

# Find all result files with the specified filter
result_files <- list.files(results_dir, pattern=paste0("results_[^_]+_", sample_filter, "\\.rda"), full.names=TRUE)
cat("Found", length(result_files), "result files to analyze with filter:", sample_filter, "\n")

# Initialize data frame for detailed results
all_detailed_results <- data.frame()

# Process each deconvolution method
for (result_file in result_files) {
  # Extract method name, accounting for the filter suffix
  method_name <- basename(result_file)
  method_name <- gsub("results_", "", method_name)
  method_name <- gsub(paste0("_", sample_filter, "\\.rda$"), "", method_name)
  
  cat("Processing method:", method_name, "\n")
  
  # Load deconvolution results
  load(result_file)
  
  
  # Extract proportion matrix
  if (!exists("deconvolutionResult") || is.null(deconvolutionResult[[method_name]]) || 
      is.null(deconvolutionResult[[method_name]]$P)) {
    cat("Skipping", method_name, "- No proportion matrix found\n")
    next
  }
  
  proportions <- deconvolutionResult[[method_name]]$P
  
  # Process each sample's predictions and match with ground truth
  sample_names <- rownames(proportions)
  matched_gt_names <- rownames(per_sample_gt_proportions)
  
  # Find matching pairs
  matched_pairs <- list()
  
  for (sample in sample_names) {
    # Check if there's a direct match
    if (sample %in% matched_gt_names) {
      matched_pairs[[sample]] <- sample
    } 
    # Check if there's a match with "-GT" suffix
    else if (paste0(sample, "-GT") %in% matched_gt_names) {
      matched_pairs[[sample]] <- paste0(sample, "-GT")
    }
  }
  
  cat("Found", length(matched_pairs), "matching sample pairs for", method_name, "\n")
  
  if (length(matched_pairs) == 0) {
    cat("Skipping", method_name, "- No matching samples found\n")
    next
  }
  
  # Initialize metrics storage
  metrics_by_sample <- data.frame()
  
  for (pred_sample in names(matched_pairs)) {
    gt_sample <- matched_pairs[[pred_sample]]
    
    # Get predicted and ground truth proportions
    pred_props <- proportions[pred_sample, , drop=FALSE]
    gt_props <- per_sample_gt_proportions[gt_sample, , drop=FALSE]
    
    # Get all unique cell types
    all_cell_types <- unique(c(colnames(pred_props), colnames(gt_props)))
    
    # Create full vectors with zeros for missing cell types
    pred_vector <- rep(0, length(all_cell_types))
    names(pred_vector) <- all_cell_types
    existing_pred_cells <- intersect(names(pred_vector), colnames(pred_props))
    pred_vector[existing_pred_cells] <- as.numeric(pred_props[1, existing_pred_cells])
    
    gt_vector <- rep(0, length(all_cell_types))
    names(gt_vector) <- all_cell_types
    existing_gt_cells <- intersect(names(gt_vector), colnames(gt_props))
    gt_vector[existing_gt_cells] <- as.numeric(gt_props[1, existing_gt_cells])
    
    # Compute metrics
    pearson_cor <- cor(pred_vector, gt_vector, method = "pearson")
    spearman_cor <- cor(pred_vector, gt_vector, method = "spearman")
    mae <- mean(abs(pred_vector - gt_vector))
    rmse <- sqrt(mean((pred_vector - gt_vector)^2))
    r_squared <- 1 - (sum((gt_vector - pred_vector)^2) / sum((gt_vector - mean(gt_vector))^2))
    
    # Calculate Jensen-Shannon Divergence
    # Note: JSD requires positive values that sum to 1
    jsd <- function(p, q) {
      # Ensure no zeros that would cause log(0)
      p <- p + 1e-10
      q <- q + 1e-10
      # Normalize
      p <- p / sum(p)
      q <- q / sum(q)
      m <- 0.5 * (p + q)
      return(0.5 * sum(p * log(p/m)) + 0.5 * sum(q * log(q/m)))
    }
    
    js_divergence <- jsd(pred_vector, gt_vector)
    
    # Record metrics for this sample
    sample_metrics <- data.frame(
      Method = method_name,
      Sample = pred_sample,
      GT_Sample = gt_sample,
      PearsonCorr = pearson_cor,
      SpearmanCorr = spearman_cor,
      MAE = mae,
      RMSE = rmse,
      R2 = r_squared,
      JSD = js_divergence,
      CellTypeCount = length(all_cell_types),
      GT_CellTypeCount = sum(gt_vector > 0),
      Pred_CellTypeCount = sum(pred_vector > 0)
    )
    
    metrics_by_sample <- rbind(metrics_by_sample, sample_metrics)
  }
  
  # Add to overall results
  all_detailed_results <- rbind(all_detailed_results, metrics_by_sample)
}

if (nrow(all_detailed_results) > 0) {
  # Save detailed results with sample filter in filename
  detailed_output_file <- file.path(output_dir, paste0(dataset_prefix, "_detailed_metrics_", sample_filter, ".csv"))
  write.csv(all_detailed_results, detailed_output_file, row.names = FALSE)
  
  cat("Detailed metrics saved to:", detailed_output_file, "\n")
  cat("Number of methods analyzed:", length(unique(all_detailed_results$Method)), "\n")
} else {
  cat("No valid results found to save\n")
}