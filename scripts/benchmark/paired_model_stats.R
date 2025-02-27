#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
  stop("5 arguments must be supplied: <R_library_path> <per_sample_ground_truth_path> <general_ground_truth_path> <results_path> <output_dir>", call.=FALSE)
}

####### Parameters
path_Rlibrary <- args[1]
per_sample_gt_path <- args[2]
general_gt_path <- args[3]
results_path <- args[4]
output_dir <- args[5]

# Libraries
.libPaths(path_Rlibrary, FALSE)
library(Metrics)  # For MAE, RMSE, R²
library(stats)    # For correlation

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load deconvolution results
load(results_path)

# Extract method name from results file path
method_name <- basename(results_path)
method_name <- gsub("results_", "", method_name)
method_name <- gsub("\\.rda$", "", method_name)
parts <- strsplit(method_name, "_")[[1]]
method_name <- parts[1]
data_name <- paste(parts[-1], collapse="_")

# Extract the proportions matrix
if (!exists("deconvolutionResult")) {
  stop("Expected 'deconvolutionResult' object not found in results file")
}

proportions <- deconvolutionResult[[method_name]]$P
if (is.null(proportions)) {
  stop(paste("No proportion matrix found for method:", method_name))
}

# Load per-sample ground truth
load(per_sample_gt_path)
per_sample_gt <- groundTruth$P
print(paste("Per-sample ground truth contains", nrow(per_sample_gt), "samples"))

# Load general ground truth for samples without pairs
load(general_gt_path)
general_gt <- groundTruth$P
print(paste("General ground truth contains", nrow(general_gt), "rows"))

# Function to find the best matching sample name
find_matching_sample <- function(pred_sample, gt_samples) {
  # Try exact match
  if (pred_sample %in% gt_samples) return(pred_sample)
  
  # Try removing common prefixes/suffixes from pred_sample
  pred_clean <- gsub("^(sample_|s_|patient_|p_)", "", pred_sample, ignore.case=TRUE)
  for (gt_sample in gt_samples) {
    gt_clean <- gsub("^(sample_|s_|patient_|p_)", "", gt_sample, ignore.case=TRUE)
    if (pred_clean == gt_clean) return(gt_sample)
  }
  
  # Try matching by sample number (assuming format like "01A", "01B")
  if (grepl("^[0-9]+[A-Za-z]$", pred_sample)) {
    sample_num <- gsub("[A-Za-z]$", "", pred_sample)
    matching_samples <- grep(paste0("^", sample_num), gt_samples)
    if (length(matching_samples) > 0) {
      return(gt_samples[matching_samples[1]])
    }
  }
  
  # No match found
  return(NULL)
}

# Initialize data structures for results
paired_metrics <- list()
unpaired_metrics <- list()
sample_mapping <- data.frame(
  PredictionSample = character(),
  GroundTruthSample = character(),
  IsPaired = logical(),
  stringsAsFactors = FALSE
)

# Process each prediction sample
for (i in 1:nrow(proportions)) {
  pred_sample <- rownames(proportions)[i]
  predicted <- as.numeric(proportions[i, ])
  
  # Find matching GT sample
  matched_sample <- find_matching_sample(pred_sample, rownames(per_sample_gt))
  
  if (!is.null(matched_sample)) {
    # We have a paired ground truth
    actual <- as.numeric(per_sample_gt[matched_sample, names(predicted)])
    
    # Record the mapping
    sample_mapping <- rbind(sample_mapping, data.frame(
      PredictionSample = pred_sample,
      GroundTruthSample = matched_sample,
      IsPaired = TRUE,
      stringsAsFactors = FALSE
    ))
    
    # Calculate metrics
    pearson_corr <- cor(predicted, actual, method = "pearson")
    spearman_corr <- cor(predicted, actual, method = "spearman")
    mae_value <- mae(actual, predicted)
    rmse_value <- rmse(actual, predicted)
    r_squared <- 1 - (sum((actual - predicted)^2) / sum((actual - mean(actual))^2))
    
    paired_metrics[[pred_sample]] <- c(
      Sample = pred_sample,
      GTSample = matched_sample,
      Pearson = pearson_corr,
      Spearman = spearman_corr,
      MAE = mae_value,
      RMSE = rmse_value,
      R2 = r_squared,
      Paired = TRUE
    )
  } else {
    # No pair found, use general ground truth
    if (ncol(general_gt) != length(predicted)) {
      # Make sure column order matches
      if (all(colnames(proportions) %in% colnames(general_gt))) {
        actual <- as.numeric(general_gt[1, colnames(proportions)])
      } else {
        warning(paste("Cannot match columns between prediction and general ground truth for sample", 
                      pred_sample))
        next
      }
    } else {
      actual <- as.numeric(general_gt)
    }
    
    # Record the mapping
    sample_mapping <- rbind(sample_mapping, data.frame(
      PredictionSample = pred_sample,
      GroundTruthSample = "general",
      IsPaired = FALSE,
      stringsAsFactors = FALSE
    ))
    
    # Calculate metrics
    pearson_corr <- cor(predicted, actual, method = "pearson")
    spearman_corr <- cor(predicted, actual, method = "spearman")
    mae_value <- mae(actual, predicted)
    rmse_value <- rmse(actual, predicted)
    r_squared <- 1 - (sum((actual - predicted)^2) / sum((actual - mean(actual))^2))
    
    unpaired_metrics[[pred_sample]] <- c(
      Sample = pred_sample,
      GTSample = "general",
      Pearson = pearson_corr,
      Spearman = spearman_corr,
      MAE = mae_value,
      RMSE = rmse_value,
      R2 = r_squared,
      Paired = FALSE
    )
  }
}

# Combine paired and unpaired metrics
all_metrics <- c(paired_metrics, unpaired_metrics)
metrics_df <- as.data.frame(do.call(rbind, all_metrics))

# Convert metrics to numeric
numeric_columns <- c("Pearson", "Spearman", "MAE", "RMSE", "R2")
metrics_df[, numeric_columns] <- lapply(metrics_df[, numeric_columns], as.numeric)
metrics_df$Paired <- as.logical(metrics_df$Paired)

# Calculate composite score as in the original script
composite_score <- function(metrics_df) {
  # Weight parameters
  w_pearson <- 1.5
  w_spearman <- 1.5
  w_mae <- 1.5
  w_rmse <- 1.5
  w_r2 <- 0.25
  
  # Score calculation
  pearson <- abs(metrics_df$Pearson)
  spearman <- abs(metrics_df$Spearman)
  r2 <- metrics_df$R2
  
  score <- w_pearson * pearson + 
           w_spearman * spearman -
           w_mae * metrics_df$MAE -
           w_rmse * metrics_df$RMSE +
           (w_r2 * r2 * abs(r2))
  
  return(score)
}

# Add composite score to metrics dataframe
metrics_df$CompositeScore <- composite_score(metrics_df)
numeric_columns <- c(numeric_columns, "CompositeScore")

# Create separate summaries for paired and unpaired samples
paired_summary <- if (sum(metrics_df$Paired) > 0) {
  data.frame(
    Metric = numeric_columns,
    Mean = colMeans(metrics_df[metrics_df$Paired, numeric_columns], na.rm = TRUE),
    SD = apply(metrics_df[metrics_df$Paired, numeric_columns], 2, sd, na.rm = TRUE),
    Type = "Paired"
  )
} else {
  data.frame(
    Metric = numeric_columns,
    Mean = NA,
    SD = NA,
    Type = "Paired"
  )
}

unpaired_summary <- if (sum(!metrics_df$Paired) > 0) {
  data.frame(
    Metric = numeric_columns,
    Mean = colMeans(metrics_df[!metrics_df$Paired, numeric_columns], na.rm = TRUE),
    SD = apply(metrics_df[!metrics_df$Paired, numeric_columns], 2, sd, na.rm = TRUE),
    Type = "Unpaired"
  )
} else {
  data.frame(
    Metric = numeric_columns,
    Mean = NA,
    SD = NA,
    Type = "Unpaired"
  )
}

overall_summary <- data.frame(
  Metric = numeric_columns,
  Mean = colMeans(metrics_df[, numeric_columns], na.rm = TRUE),
  SD = apply(metrics_df[, numeric_columns], 2, sd, na.rm = TRUE),
  Type = "Overall"
)

summary_stats <- rbind(paired_summary, unpaired_summary, overall_summary)

# Save results
metrics_file <- file.path(output_dir, paste0(method_name, "_", data_name, "_paired_benchmarking.csv"))
write.csv(metrics_df, metrics_file, row.names = FALSE)
print(paste("Paired benchmarking results saved to:", metrics_file))

mapping_file <- file.path(output_dir, paste0(method_name, "_", data_name, "_sample_mapping.csv"))
write.csv(sample_mapping, mapping_file, row.names = FALSE)
print(paste("Sample mapping saved to:", mapping_file))

summary_file <- file.path(output_dir, paste0(method_name, "_", data_name, "_paired_summary_stats.csv"))
write.csv(summary_stats, summary_file, row.names = FALSE)
print(paste("Paired summary statistics saved to:", summary_file))