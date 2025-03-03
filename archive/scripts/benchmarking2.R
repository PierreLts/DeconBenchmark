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
library(ggplot2)
library(reshape2)
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

# Visualization (unchanged)
gt_df <- as.data.frame(gt_proportions)
gt_df$Sample <- "Ground Truth"  
spacer_df <- gt_df
spacer_df$Sample <- " "  # Spacer for visual clarity

# Convert to long format
gt_long <- melt(gt_df, id.vars="Sample", variable.name="CellType", value.name="Proportion")
spacer_long <- melt(spacer_df, id.vars="Sample", variable.name="CellType", value.name="Proportion")
spacer_long$Proportion <- 0 

# Convert predicted proportions to long format
prop_df <- as.data.frame(proportions)
prop_df$Sample <- rownames(proportions)
prop_long <- melt(prop_df, id.vars="Sample", variable.name="CellType", value.name="Proportion")

# Combine predicted, spacer, and ground truth data
combined_data <- rbind(prop_long, spacer_long, gt_long)

# Convert Sample to factor
combined_data$Sample <- factor(combined_data$Sample, levels = c(unique(prop_long$Sample), " ", "Ground Truth"))

# Generate color palette
num_cell_types <- length(unique(combined_data$CellType))
custom_colors <- colorRampPalette(c("#E57373", "#FFB74D", "#FFD54F", "#AED581", "#4FC3F7", "#64B5F6", "#BA68C8", "#F06292", "#FF80AB", "#7986CB"))(num_cell_types)
names(custom_colors) <- unique(combined_data$CellType)

# Create stacked bar plot
p <- ggplot(combined_data, aes(x=Sample, y=Proportion, fill=CellType)) +
  geom_bar(stat="identity", position="stack") +
  theme_minimal() +
  labs(title=paste(method_name, "vs Ground Truth Cell Type Proportions"),
       x="Samples", y="Proportions") +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        panel.grid.major=element_line(color="lightgray"),
        panel.grid.minor=element_line(color="lightgray"),
        legend.position="right",
        legend.key.size = unit(0.5, "cm"),
        plot.title=element_text(hjust=0.5)) +
  scale_fill_manual(values=custom_colors, name="Cell Type")

# Save the plot
output_filename <- file.path(output_dir, paste0(method_name, "_", data_name, "_vs_ground_truth.pdf"))
ggsave(output_filename, p, width=12, height=6)
print(paste("Plot saved to:", output_filename))