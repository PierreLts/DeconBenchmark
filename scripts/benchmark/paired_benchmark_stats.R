#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop(paste("4 arguments must be supplied instead of", length(args)), call. = FALSE)
}

####### Parameter of script
path_Rlibrary <- args[1] #IMPORTANT
results_dir <- args[2]    # Directory containing deconvolution results
ground_truth_path <- args[3]  # Path to per-sample ground truth file
output_dir <- args[4]     # Output directory for results

# Libraries
.libPaths(path_Rlibrary, FALSE) #IMPORTANT
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(Metrics)  # For common error metrics
library(scales)   # For rescaling data
print("Check Libraries")

# Load per-sample ground truth
load(ground_truth_path)
per_sample_gt_proportions <- groundTruth$P

# Find all result files
result_files <- list.files(results_dir, pattern="results_.*\\.rda", full.names=TRUE)
cat("Found", length(result_files), "result files to analyze\n")

# Initialize data frames for summary results
all_method_summaries <- data.frame()
all_detailed_results <- data.frame()

# Process each deconvolution method
for (result_file in result_files) {
  method_name <- basename(result_file)
  method_name <- gsub("results_", "", method_name)
  method_name <- gsub("\\.rda$", "", method_name)
  
  # Extract main method name from potentially composite name
  parts <- strsplit(method_name, "_")[[1]]
  method_name <- parts[1]
  data_name <- paste(parts[-1], collapse="_")
  
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

# Save results
if (nrow(all_method_summaries) > 0) {
  # Calculate composite score
  # First normalize all metrics to 0-1 scale
  normalize <- function(x, higher_is_better = TRUE) {
    if (length(unique(x)) == 1) return(rep(0.5, length(x)))  # Handle constant values
    if (higher_is_better) {
      return(scales::rescale(x, to = c(0, 1)))
    } else {
      return(scales::rescale(-x, to = c(0, 1)))  # Flip for lower-is-better metrics
    }
  }
  
  all_method_summaries <- all_method_summaries %>%
    mutate(
      norm_pearson = normalize(Mean_PearsonCorr, TRUE),
      norm_spearman = normalize(Mean_SpearmanCorr, TRUE),
      norm_mae = normalize(Mean_MAE, FALSE),
      norm_rmse = normalize(Mean_RMSE, FALSE),
      norm_r2 = normalize(Mean_R2, TRUE),
      CompositeScore = (norm_pearson + norm_spearman + norm_mae + norm_rmse + norm_r2) / 5
    )
  
  # Rank methods by composite score
  ranked_methods <- all_method_summaries %>%
    arrange(desc(CompositeScore))
  
  # Save detailed results
  detailed_output_file <- file.path(output_dir, "paired_sample_detailed_metrics.csv")
  write.csv(all_detailed_results, detailed_output_file, row.names = FALSE)
  
  # Save summary results
  summary_output_file <- file.path(output_dir, "paired_sample_method_summary.csv")
  write.csv(ranked_methods, summary_output_file, row.names = FALSE)
  
  # Generate comparison plot with all metrics
  if (nrow(ranked_methods) > 1) {
    # Select top methods (limit to make plot readable)
    top_n_methods <- min(13, nrow(ranked_methods))
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
    
    # Set up color palette for metrics - match the image
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
    
    # Create plot
    p <- ggplot(plot_data_long, aes(x = RankedName, y = MetricValue, fill = MetricType)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
      geom_errorbar(aes(ymin = MetricValue - ErrorValue, ymax = MetricValue + ErrorValue), 
                   position = position_dodge(width = 0.8), width = 0.25) +
      labs(
        title = "Performance Metrics Comparison Across Deconvolution Methods",
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
      scale_y_continuous(limits = c(-15, 5), breaks = seq(-15, 5, 5))
    
    # Save plot
    comparison_plot_file <- file.path(output_dir, "method_metrics_comparison.pdf")
    ggsave(comparison_plot_file, p, width = 22, height = 8)
    
    # Also save as PNG for easier viewing
    comparison_plot_png <- file.path(output_dir, "method_metrics_comparison.png")
    ggsave(comparison_plot_png, p, width = 18, height = 8, dpi = 300)
    
    cat("Comparison plot saved to:", comparison_plot_file, "and", comparison_plot_png, "\n")
  }
  
  cat("Results saved to:", output_dir, "\n")
  cat("Top performing methods by Composite Score:\n")
  print(head(ranked_methods[, c("Method", "CompositeScore", "Mean_PearsonCorr", "Mean_MAE", "SampleCount")], 5))
} else {
  cat("No valid results found to save\n")
}