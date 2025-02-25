#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop(paste("4 arguments must be supplied instead of", length(args)), call. = FALSE)
}

####### Script Parameters (ORDER IS IMPORTANT)
path_Rlibrary <- args[1] # IMPORTANT
ground_truth_path <- args[2]
results_path <- args[3]
output_dir <- args[4]

# Libraries
.libPaths(path_Rlibrary, FALSE) # IMPORTANT
library(Metrics)  # For MAE, RMSE, R²
library(stats)    # For correlation

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

# Load ground truth
load(ground_truth_path)
gt_proportions <- groundTruth$P

# Ensure ground truth has the same cell types as predicted values
if (!all(colnames(proportions) %in% colnames(gt_proportions))) {
  stop("Mismatch: Ground truth and predicted proportions do not have the same cell types.")
}

# **Fix dimension mismatch**
if (nrow(gt_proportions) == 1) {
  gt_proportions <- gt_proportions[rep(1, nrow(proportions)), ]  # Replicate for each sample
} else if (nrow(gt_proportions) != nrow(proportions)) {
  stop("Mismatch: Predicted and ground truth proportions have different number of samples.")
}

# Compute **benchmarking metrics per sample**
metrics_list <- list()
for (i in 1:nrow(proportions)) {
  predicted <- as.numeric(proportions[i, ])
  actual <- as.numeric(gt_proportions[i, ])
  
  pearson_corr <- cor(predicted, actual, method = "pearson")
  spearman_corr <- cor(predicted, actual, method = "spearman")
  mae_value <- mae(actual, predicted)
  rmse_value <- rmse(actual, predicted)
  r_squared <- 1 - (sum((actual - predicted)^2) / sum((actual - mean(actual))^2))

  metrics_list[[i]] <- c(
    Sample = rownames(proportions)[i],
    Pearson = pearson_corr,
    Spearman = spearman_corr,
    MAE = mae_value,
    RMSE = rmse_value,
    R2 = r_squared
  )
}

# Convert metrics_df to a proper numeric format
metrics_df <- as.data.frame(do.call(rbind, metrics_list))
metrics_df$Sample <- as.character(metrics_df$Sample)  # Ensure Sample is a string

# Convert numeric columns explicitly
numeric_columns <- c("Pearson", "Spearman", "MAE", "RMSE", "R2")
metrics_df[, numeric_columns] <- lapply(metrics_df[, numeric_columns], as.numeric)

# Calculate composite score
composite_score <- function(metrics_df) {
  # Weight parameters
  w_pearson <- 1.5
  w_spearman <- 1.5
  w_mae <- 1.5
  w_rmse <- 1.5
  w_r2 <- 0.25
  
  # Normalize MAE and RMSE (lower is better)
  pearson <- abs(metrics_df$Pearson)
  spearman <- abs(metrics_df$Spearman)
  r2 <- metrics_df$R2
  
  # Calculate composite score
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

# Summary statistics: Mean ± SD
summary_stats <- data.frame(
  Metric = numeric_columns,
  Mean = colMeans(metrics_df[, numeric_columns], na.rm = TRUE),
  SD = apply(metrics_df[, numeric_columns], 2, sd, na.rm = TRUE)
)

# Save metrics to file
metrics_file <- file.path(output_dir, paste0(method_name, "_", data_name, "_per_sample_benchmarking.csv"))
write.csv(metrics_df, metrics_file, row.names = FALSE)
print(paste("Per-sample benchmarking results saved to:", metrics_file))

summary_file <- file.path(output_dir, paste0(method_name, "_", data_name, "_summary_stats.csv"))
write.csv(summary_stats, summary_file, row.names = FALSE)
print(paste("Summary statistics saved to:", summary_file))
