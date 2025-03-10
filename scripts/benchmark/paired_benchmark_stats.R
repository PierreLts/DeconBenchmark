#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop(paste("4 arguments must be supplied instead of", length(args)), call. = FALSE)
}

####### Parameter of script
path_Rlibrary <- args[1] #IMPORTANT
dataset_prefix <- args[2] # Dataset prefix (e.g., TB1)
output_base_dir <- args[3] # Base output directory for benchmark results
include_overall_gt <- as.logical(args[4]) # Whether to include overall ground truth

# Libraries
.libPaths(path_Rlibrary, FALSE) #IMPORTANT
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(Metrics)  # For common error metrics
library(scales)   # For rescaling data
print("Check Libraries")

# Setup paths
results_dir <- file.path("/work/gr-fe/lorthiois/DeconBenchmark/deconv_results", dataset_prefix)
data_dir <- file.path("/work/gr-fe/lorthiois/DeconBenchmark/generated_data", dataset_prefix)
output_dir <- file.path(output_base_dir, dataset_prefix, "benchmarks")

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load per-sample ground truth
per_sample_gt_path <- file.path(data_dir, paste0(dataset_prefix, "_GT_proportions_per_sample.rda"))
if (!file.exists(per_sample_gt_path)) {
  stop(paste("Per-sample ground truth file not found:", per_sample_gt_path))
}
load(per_sample_gt_path)
per_sample_gt_proportions <- groundTruth$P

# Find all result files
result_files <- list.files(results_dir, pattern="results_[^_]+\\.rda", full.names=TRUE)
cat("Found", length(result_files), "result files to analyze\n")

# Initialize data frames for summary results
all_method_summaries <- data.frame()
all_detailed_results <- data.frame()

# Process each deconvolution method
for (result_file in result_files) {
  method_name <- basename(result_file)
  method_name <- gsub("results_", "", method_name)
  method_name <- gsub("\\.rda$", "", method_name)
  
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
  
  # Compute summary statistics for this method
  method_summary <- metrics_by_sample %>%
    summarize(
      Method = first(Method),
      SampleCount = n(),
      Mean_PearsonCorr = mean(PearsonCorr, na.rm=TRUE),
      SD_PearsonCorr = sd(PearsonCorr, na.rm=TRUE),
      Mean_SpearmanCorr = mean(SpearmanCorr, na.rm=TRUE),
      SD_SpearmanCorr = sd(SpearmanCorr, na.rm=TRUE),
      Mean_MAE = mean(MAE, na.rm=TRUE),
      SD_MAE = sd(MAE, na.rm=TRUE),
      Mean_RMSE = mean(RMSE, na.rm=TRUE),
      SD_RMSE = sd(RMSE, na.rm=TRUE),
      Mean_R2 = mean(R2, na.rm=TRUE),
      SD_R2 = sd(R2, na.rm=TRUE),
      Mean_JSD = mean(JSD, na.rm=TRUE),
      SD_JSD = sd(JSD, na.rm=TRUE)
    )
  
  # Add to overall results
  all_method_summaries <- rbind(all_method_summaries, method_summary)
  all_detailed_results <- rbind(all_detailed_results, metrics_by_sample)
}

# Define the composite score function
composite_score <- function(metrics_df) {
  # Weight parameters
  w_pearson <- 1.5
  w_spearman <- 1.5
  w_mae <- 1.5
  w_rmse <- 1.5
  w_r2 <- 0
  
  # Calculate composite score
  pearson <- abs(metrics_df$Pearson)
  spearman <- abs(metrics_df$Spearman)
  r2 <- metrics_df$R2
  
  #Calculate composite score
  score <- ((w_pearson * pearson) + 
           (w_spearman * spearman) -
           (w_mae * metrics_df$MAE) -
           (w_rmse * metrics_df$RMSE) +
           (w_r2 * r2))/4

  return(score)
}

if (nrow(all_method_summaries) > 0) {
  # Calculate composite score
  all_method_summaries <- all_method_summaries %>%
    mutate(
      CompositeScore = composite_score(data.frame(
        Pearson = Mean_PearsonCorr,
        Spearman = Mean_SpearmanCorr,
        MAE = Mean_MAE,
        RMSE = Mean_RMSE,
        R2 = Mean_R2
      ))
    )

  # Rank methods by composite score
  ranked_methods <- all_method_summaries %>%
    arrange(desc(CompositeScore))
  
  # Save detailed results
  detailed_output_file <- file.path(output_dir, paste0(dataset_prefix, "_detailed_metrics.csv"))
  write.csv(all_detailed_results, detailed_output_file, row.names = FALSE)
  
  # Save summary results
  summary_output_file <- file.path(output_dir, paste0(dataset_prefix, "_method_summary.csv"))
  write.csv(ranked_methods, summary_output_file, row.names = FALSE)
  
  # Generate comparison plot with all metrics
  if (nrow(ranked_methods) > 1) {
    # Select top methods (limit to make plot readable)
    top_n_methods <- min(15, nrow(ranked_methods))
    top_methods <- ranked_methods[1:top_n_methods, ]
    
    # Add ranking numbers to method names
    top_methods$RankedName <- paste0(1:nrow(top_methods), ". ", top_methods$Method)
    
    # Prepare data for plotting
    plot_data <- top_methods %>%
      select(Method, RankedName,
             Mean_PearsonCorr, SD_PearsonCorr, 
             Mean_SpearmanCorr, SD_SpearmanCorr,
             Mean_MAE, SD_MAE,
             Mean_RMSE, SD_RMSE,
             Mean_R2, SD_R2,
             CompositeScore)
    
    # Reshape data for plotting
    plot_data_long <- data.frame(
      Method = rep(plot_data$Method, 6),
      RankedName = rep(plot_data$RankedName, 6),
      MetricValue = c(
        plot_data$Mean_PearsonCorr,
        plot_data$Mean_SpearmanCorr,
        -plot_data$Mean_MAE,       # Negate MAE for visualization (lower is better)
        -plot_data$Mean_RMSE,      # Negate RMSE for visualization (lower is better)
        plot_data$Mean_R2,
        plot_data$CompositeScore * 10 - 5  # Scale to match other metrics
      ),
      ErrorValue = c(
        plot_data$SD_PearsonCorr,
        plot_data$SD_SpearmanCorr,
        plot_data$SD_MAE,
        plot_data$SD_RMSE,
        plot_data$SD_R2,
        rep(0.5, nrow(plot_data))  # Fixed error bar for composite
      ),
      MetricType = rep(c("Pearson (+)", "Spearman (+)", "MAE (-)", "RMSE (-)", "R² (+)", "Composite (+)"), each = nrow(plot_data))
    )
    
    # Set up color palette for metrics
    metric_colors <- c(
    "Pearson (+)" = "#4CBFFF",
    "Spearman (+)" = "#AED581",
    "MAE (-)" = "#FF8CBF",
    "RMSE (-)" = "#FF9933",
    "R² (+)" = "#BA68C8",
    "Composite (+)" = "#E64CA6"
    )
    
    # Order methods by composite score (best to worst)
    plot_data_long$RankedName <- factor(plot_data_long$RankedName, 
                                       levels = plot_data$RankedName)
    
    # Find adaptive y-axis limits based on data
    max_with_error <- max(abs(plot_data_long$MetricValue) + plot_data_long$ErrorValue, na.rm = TRUE)

    # Add buffer and round up
    max_abs_value <- ceiling(max_with_error * 1.02)
    # Ensure a minimum range for visibility
    max_abs_value <- max(max_abs_value, 1.0)  # Set minimum to 1.0
    # Create sensible breaks
    break_step <- max(1, round(max_abs_value/5))
    break_sequence <- seq(-max_abs_value, max_abs_value, by=break_step)

    # Create plot
    p <- ggplot(plot_data_long, aes(x = RankedName, y = MetricValue, fill = MetricType)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.85), width = 0.8) +
    geom_errorbar(aes(ymin = MetricValue - ErrorValue, ymax = MetricValue + ErrorValue), 
                position = position_dodge(width = 0.8), width = 0) +
    labs(
        title = paste("Performance Metrics Comparison -", dataset_prefix),
        subtitle = "(+) higher is better, (-) lower is better",
        x = "",
        y = "Score Value"
    ) +
    scale_fill_manual(values = metric_colors, name = "MetricDisplay") +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(color = "lightgray"),
        panel.grid.minor.y = element_line(color = "lightgray"),
        panel.grid.major.x = element_blank()
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    scale_y_continuous(
        limits = c(-max_abs_value, max_abs_value),
        breaks = break_sequence
    )

    # Save plot as PDF
    comparison_plot_file <- file.path(output_dir, paste0(dataset_prefix, "_method_metrics_comparison.pdf"))
    ggsave(comparison_plot_file, p, width = 22, height = 12)
    
    
    cat("Comparison plots saved to:", comparison_plot_file, "and", comparison_plot_png, "\n")
  }
  
  cat("Results saved to:", output_dir, "\n")
  cat("Top performing methods by Composite Score:\n")
  print(head(ranked_methods[, c("Method", "CompositeScore", "Mean_PearsonCorr", "Mean_MAE", "SampleCount")], 5))
} else {
  cat("No valid results found to save\n")
}