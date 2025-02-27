#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
  stop("5 arguments must be supplied: <R_library_path> <mean_ground_truth> <patient_ground_truth_dir> <results_path> <output_dir>", call.=FALSE)
}

####### Script Parameters 
path_Rlibrary <- args[1]
mean_ground_truth_path <- args[2]
patient_ground_truth_dir <- args[3]
results_path <- args[4]
output_dir <- args[5]

# Libraries
.libPaths(path_Rlibrary, FALSE)
library(Metrics)  # For MAE, RMSE, R²
library(stats)    # For correlation
library(dplyr)
library(reshape2)
library(ggplot2)

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

# Load mean ground truth (for samples without patient-specific matches)
load(mean_ground_truth_path)
mean_gt_proportions <- groundTruth$P

# Load patient mapping data
patient_mapping_file <- file.path(patient_ground_truth_dir, "patient_sample_mapping.csv")
if (file.exists(patient_mapping_file)) {
  patient_mapping <- read.csv(patient_mapping_file, stringsAsFactors = FALSE)
} else {
  warning("Patient mapping file not found. Using only mean ground truth.")
  patient_mapping <- data.frame(sample = character(0), patient_id = character(0))
}

# Get list of patient ground truth files
patient_gt_files <- list.files(
  path = patient_ground_truth_dir, 
  pattern = "ground_truth_patient_.*\\.rda$",
  full.names = TRUE
)

# Load all patient ground truths into a list
patient_ground_truths <- list()
for (gt_file in patient_gt_files) {
  patient_id <- gsub(".*ground_truth_patient_(.+)\\.rda$", "\\1", gt_file)
  load(gt_file)  # Loads patientGroundTruth
  patient_ground_truths[[patient_id]] <- patientGroundTruth
}

# Define weights for patient-specific vs mean ground truth
PATIENT_MATCH_WEIGHT <- 0.8  # Higher weight for patient-specific matches
MEAN_MATCH_WEIGHT <- 0.2     # Lower weight for mean-based matches

# Initialize dataframes to hold metrics 
metrics_list <- list()
patient_matched_samples <- c()
mean_matched_samples <- c()

# Process each sample in the deconvolution results
for (i in 1:nrow(proportions)) {
  sample_id <- rownames(proportions)[i]
  predicted <- as.numeric(proportions[i, ])
  
  # Check if this sample has a patient-specific match
  patient_match <- patient_mapping[patient_mapping$sample == sample_id, ]
  
  if (nrow(patient_match) > 0) {
    # We have a patient-specific match
    patient_id <- patient_match$patient_id[1]
    
    if (patient_id %in% names(patient_ground_truths)) {
      # Get the patient-specific ground truth
      patient_gt <- patient_ground_truths[[patient_id]]
      
      # Find the exact ground truth row (could be before/after treatment)
      gt_row <- which(rownames(patient_gt$P) == sample_id)
      
      if (length(gt_row) > 0) {
        actual <- as.numeric(patient_gt$P[gt_row, ])
        match_type <- "patient_specific"
        match_weight <- PATIENT_MATCH_WEIGHT
        patient_matched_samples <- c(patient_matched_samples, sample_id)
      } else {
        # Fall back to mean if exact row not found
        actual <- as.numeric(mean_gt_proportions[1, ]) # Using first row as example
        match_type <- "mean_fallback"
        match_weight <- MEAN_MATCH_WEIGHT
        mean_matched_samples <- c(mean_matched_samples, sample_id)
      }
    } else {
      # Patient ground truth not found
      actual <- as.numeric(mean_gt_proportions[1, ])
      match_type <- "mean_fallback"
      match_weight <- MEAN_MATCH_WEIGHT
      mean_matched_samples <- c(mean_matched_samples, sample_id)
    }
  } else {
    # No patient match found, use mean ground truth
    actual <- as.numeric(mean_gt_proportions[1, ])
    match_type <- "mean_only"
    match_weight <- MEAN_MATCH_WEIGHT
    mean_matched_samples <- c(mean_matched_samples, sample_id)
  }
  
  # Ensure ground truth has the same cell types
  if (!all(colnames(proportions) %in% names(actual))) {
    # Need to reorder or subset cell types
    actual_aligned <- rep(0, length(predicted))
    names(actual_aligned) <- colnames(proportions)
    
    for (ct in names(actual)) {
      if (ct %in% colnames(proportions)) {
        actual_aligned[ct] <- actual[ct]
      }
    }
    actual <- actual_aligned
  }
  
  # Calculate metrics
  pearson_corr <- cor(predicted, actual, method = "pearson")
  spearman_corr <- cor(predicted, actual, method = "spearman")
  mae_value <- mae(actual, predicted)
  rmse_value <- rmse(actual, predicted)
  r_squared <- 1 - (sum((actual - predicted)^2) / sum((actual - mean(actual))^2))
  
  metrics_list[[i]] <- c(
    Sample = sample_id,
    Pearson = pearson_corr,
    Spearman = spearman_corr,
    MAE = mae_value,
    RMSE = rmse_value,
    R2 = r_squared,
    MatchType = match_type,
    Weight = match_weight
  )
}

# Convert metrics to a data frame
metrics_df <- as.data.frame(do.call(rbind, metrics_list))
metrics_df$Sample <- as.character(metrics_df$Sample)

# Convert numeric columns
numeric_columns <- c("Pearson", "Spearman", "MAE", "RMSE", "R2", "Weight")
metrics_df[, numeric_columns] <- lapply(metrics_df[, numeric_columns], as.numeric)

# Calculate composite score with weighting
composite_score <- function(metrics_df) {
  # Weight parameters
  w_pearson <- 1.5
  w_spearman <- 1.5
  w_mae <- 1.5
  w_rmse <- 1.5
  w_r2 <- 0.25
  
  # Calculate basic score
  pearson <- abs(metrics_df$Pearson)
  spearman <- abs(metrics_df$Spearman)
  r2 <- metrics_df$R2
  
  # Calculate composite score
  score <- w_pearson * pearson + 
           w_spearman * spearman -
           w_mae * metrics_df$MAE -
           w_rmse * metrics_df$RMSE +
           (w_r2 * r2 * abs(r2))
  
  # Apply sample match weight
  weighted_score <- score * metrics_df$Weight
  
  return(weighted_score)
}

# Add weighted composite score to metrics
metrics_df$WeightedScore <- composite_score(metrics_df)

# Print summary of matching
print(paste("Total samples:", nrow(metrics_df)))
print(paste("Patient-matched samples:", length(patient_matched_samples), 
            "(", round(100*length(patient_matched_samples)/nrow(metrics_df), 1), "%)"))
print(paste("Mean-matched samples:", length(mean_matched_samples),
            "(", round(100*length(mean_matched_samples)/nrow(metrics_df), 1), "%)"))

# Calculate weighted summary statistics
summary_stats <- metrics_df %>%
  group_by(MatchType) %>%
  summarize(
    SampleCount = n(),
    WeightValue = first(Weight),
    Pearson_Mean = mean(Pearson, na.rm = TRUE),
    Pearson_SD = sd(Pearson, na.rm = TRUE),
    Spearman_Mean = mean(Spearman, na.rm = TRUE),
    Spearman_SD = sd(Spearman, na.rm = TRUE),
    MAE_Mean = mean(MAE, na.rm = TRUE),
    MAE_SD = sd(MAE, na.rm = TRUE),
    RMSE_Mean = mean(RMSE, na.rm = TRUE),
    RMSE_SD = sd(RMSE, na.rm = TRUE),
    R2_Mean = mean(R2, na.rm = TRUE),
    R2_SD = sd(R2, na.rm = TRUE),
    WeightedScore_Mean = mean(WeightedScore, na.rm = TRUE),
    WeightedScore_SD = sd(WeightedScore, na.rm = TRUE)
  )

# Calculate overall weighted means
overall_stats <- data.frame(
  MatchType = "overall_weighted",
  SampleCount = nrow(metrics_df),
  WeightValue = NA,
  Pearson_Mean = weighted.mean(metrics_df$Pearson, metrics_df$Weight, na.rm = TRUE),
  Pearson_SD = NA,
  Spearman_Mean = weighted.mean(metrics_df$Spearman, metrics_df$Weight, na.rm = TRUE),
  Spearman_SD = NA,
  MAE_Mean = weighted.mean(metrics_df$MAE, metrics_df$Weight, na.rm = TRUE),
  MAE_SD = NA,
  RMSE_Mean = weighted.mean(metrics_df$RMSE, metrics_df$Weight, na.rm = TRUE),
  RMSE_SD = NA,
  R2_Mean = weighted.mean(metrics_df$R2, metrics_df$Weight, na.rm = TRUE),
  R2_SD = NA,
  WeightedScore_Mean = weighted.mean(metrics_df$WeightedScore, metrics_df$Weight, na.rm = TRUE),
  WeightedScore_SD = NA
)

# Combine summaries
summary_stats <- rbind(summary_stats, overall_stats)

# Save metrics to files
metrics_file <- file.path(output_dir, paste0(method_name, "_", data_name, "_patient_benchmarking.csv"))
write.csv(metrics_df, metrics_file, row.names = FALSE)
print(paste("Patient benchmarking results saved to:", metrics_file))

summary_file <- file.path(output_dir, paste0(method_name, "_", data_name, "_patient_summary.csv"))
write.csv(summary_stats, summary_file, row.names = FALSE)
print(paste("Patient summary statistics saved to:", summary_file))

# Create visualizations
# Boxplot of metrics by match type
plot_data <- melt(metrics_df, 
                 id.vars = c("Sample", "MatchType", "Weight"),
                 measure.vars = c("Pearson", "Spearman", "MAE", "RMSE", "R2", "WeightedScore"))

p <- ggplot(plot_data, aes(x = MatchType, y = value, fill = MatchType)) +
  geom_boxplot() +
  facet_wrap(~ variable, scales = "free_y") +
  theme_minimal() +
  labs(title = paste(method_name, "Metrics by Match Type"),
       x = "Match Type", 
       y = "Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save plot
plot_file <- file.path(output_dir, paste0(method_name, "_", data_name, "_patient_metrics.pdf"))
ggsave(plot_file, p, width = 10, height = 8)
print(paste("Metrics plot saved to:", plot_file))

# Create comparison bar chart for overall vs mean-only performance
if ("patient_specific" %in% summary_stats$MatchType && "mean_only" %in% summary_stats$MatchType) {
  compare_data <- summary_stats %>%
    filter(MatchType %in% c("patient_specific", "mean_only", "overall_weighted")) %>%
    select(MatchType, Pearson_Mean, Spearman_Mean, MAE_Mean, RMSE_Mean, R2_Mean) %>%
    melt(id.vars = "MatchType")
  
  p2 <- ggplot(compare_data, aes(x = variable, y = value, fill = MatchType)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    labs(title = paste(method_name, "Performance Comparison"),
         x = "Metric", 
         y = "Value") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save comparison plot
  compare_file <- file.path(output_dir, paste0(method_name, "_", data_name, "_comparison.pdf"))
  ggsave(compare_file, p2, width = 8, height = 6)
  print(paste("Comparison plot saved to:", compare_file))
}