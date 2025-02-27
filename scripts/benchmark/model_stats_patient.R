#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop("4 arguments must be supplied: <R_library_path> <patient_ground_truth> <results_path> <output_dir>", call.=FALSE)
}

####### Script Parameters 
path_Rlibrary <- args[1]
patient_ground_truth_path <- args[2]
results_path <- args[3]
output_dir <- args[4]

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

# Load patient ground truth data
load(patient_ground_truth_path)  # Loads 'groundTruth' object with P matrix and patient_mapping

# Define weights
PATIENT_MATCH_WEIGHT <- 0.8  # Higher weight for patient-specific matches
MEAN_MATCH_WEIGHT <- 0.2     # Lower weight for mean-based matches

# Calculate overall mean ground truth for samples without exact matches
mean_gt_proportions <- colMeans(groundTruth$P)

# Initialize dataframes to hold metrics 
metrics_list <- list()

# Process each sample in the deconvolution results
for (i in 1:nrow(proportions)) {
  sample_id <- rownames(proportions)[i]
  predicted <- as.numeric(proportions[i, ])
  
  # Check if this sample has a match in the ground truth
  if (sample_id %in% rownames(groundTruth$P)) {
    # We have an exact match
    actual <- as.numeric(groundTruth$P[sample_id, ])
    match_type <- "exact_match"
    match_weight <- PATIENT_MATCH_WEIGHT
    
    # Get patient ID for this sample
    patient_id <- groundTruth$patient_mapping$patient_id[groundTruth$patient_mapping$sample == sample_id]
  } else {
    # No exact match, try to use overall mean
    actual <- mean_gt_proportions
    match_type <- "mean_based"
    match_weight <- MEAN_MATCH_WEIGHT
    patient_id <- "unknown"
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
    PatientID = patient_id,
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
metrics_df$PatientID <- as.character(metrics_df$PatientID)
metrics_df$MatchType <- as.character(metrics_df$MatchType)

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
print(paste("Exact matched samples:", sum(metrics_df$MatchType == "exact_match"), 
            "(", round(100*sum(metrics_df$MatchType == "exact_match")/nrow(metrics_df), 1), "%)"))
print(paste("Mean-based samples:", sum(metrics_df$MatchType == "mean_based"),
            "(", round(100*sum(metrics_df$MatchType == "mean_based")/nrow(metrics_df), 1), "%)"))

# Calculate per-match-type summary statistics
match_summary <- metrics_df %>%
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

# Calculate per-patient summary statistics
patient_summary <- metrics_df %>%
  filter(PatientID != "unknown") %>%
  group_by(PatientID) %>%
  summarize(
    SampleCount = n(),
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

# Calculate overall weighted statistics
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
match_summary <- rbind(match_summary, overall_stats)

# Save metrics to files
metrics_file <- file.path(output_dir, paste0(method_name, "_", data_name, "_patient_benchmarking.csv"))
write.csv(metrics_df, metrics_file, row.names = FALSE)
print(paste("Patient benchmarking results saved to:", metrics_file))

match_summary_file <- file.path(output_dir, paste0(method_name, "_", data_name, "_match_summary.csv"))
write.csv(match_summary, match_summary_file, row.names = FALSE)
print(paste("Match type summary saved to:", match_summary_file))

patient_summary_file <- file.path(output_dir, paste0(method_name, "_", data_name, "_patient_summary.csv"))
write.csv(patient_summary, patient_summary_file, row.names = FALSE)
print(paste("Patient summary saved to:", patient_summary_file))

# Create visualizations
# 1. Box plot of metrics by match type
plot_data <- melt(metrics_df, 
                 id.vars = c("Sample", "PatientID", "MatchType", "Weight"),
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
boxplot_file <- file.path(output_dir, paste0(method_name, "_", data_name, "_match_metrics.pdf"))
ggsave(boxplot_file, p, width = 10, height = 8)
print(paste("Match metrics plot saved to:", boxplot_file))

# 2. Patient variation - show how each patient performs
if (nrow(patient_summary) > 1) {
  # Prepare data for plotting
  patient_plot_data <- melt(patient_summary, 
                          id.vars = c("PatientID", "SampleCount"),
                          measure.vars = c("Pearson_Mean", "Spearman_Mean", "MAE_Mean", "RMSE_Mean", "R2_Mean"))
  
  # Rename variables for cleaner labels
  patient_plot_data$variable <- gsub("_Mean$", "", patient_plot_data$variable)
  
  # Create the plot
  p2 <- ggplot(patient_plot_data, aes(x = PatientID, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~ variable, scales = "free_y") +
    theme_minimal() +
    labs(title = paste(method_name, "Performance by Patient"),
         x = "Patient ID", 
         y = "Mean Value") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save patient plot
  patient_file <- file.path(output_dir, paste0(method_name, "_", data_name, "_patient_variation.pdf"))
  ggsave(patient_file, p2, width = 12, height = 8)
  print(paste("Patient variation plot saved to:", patient_file))
}

print("Analysis complete!")