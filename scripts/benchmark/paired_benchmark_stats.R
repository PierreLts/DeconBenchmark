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
      Mean_JSD = mean(JSD, na.rm=TRUE),
      SD_JSD = sd(JSD, na.rm=TRUE)
    )
  
  # Add to overall results
  all_method_summaries <- rbind(all_method_summaries, method_summary)
  all_detailed_results <- rbind(all_detailed_results, metrics_by_sample)
}

# Save results
if (nrow(all_method_summaries) > 0) {
  # Rank methods
  ranked_methods <- all_method_summaries %>%
    arrange(desc(Mean_PearsonCorr)) %>%
    mutate(Rank_Pearson = row_number())
  
  # Save detailed results
  detailed_output_file <- file.path(output_dir, "paired_sample_detailed_metrics.csv")
  write.csv(all_detailed_results, detailed_output_file, row.names = FALSE)
  
  # Save summary results
  summary_output_file <- file.path(output_dir, "paired_sample_method_summary.csv")
  write.csv(ranked_methods, summary_output_file, row.names = FALSE)
  
  # Generate plots
  # 1. Method comparison plot
  if (nrow(ranked_methods) > 1) {
    p1 <- ggplot(ranked_methods, aes(x=reorder(Method, Mean_PearsonCorr))) +
      geom_bar(aes(y=Mean_PearsonCorr), stat="identity", fill="steelblue") +
      geom_errorbar(aes(ymin=Mean_PearsonCorr-SD_PearsonCorr, 
                        ymax=Mean_PearsonCorr+SD_PearsonCorr), width=.2) +
      labs(title="Method Comparison: Mean Pearson Correlation",
           x="Method", y="Pearson Correlation") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle=45, hjust=1))
    
    ggsave(file.path(output_dir, "method_comparison_pearson.pdf"), p1, width=10, height=6)
    
    # 2. Error metrics comparison
    p2 <- ggplot(ranked_methods, aes(x=reorder(Method, -Mean_MAE))) +
      geom_bar(aes(y=Mean_MAE), stat="identity", fill="salmon") +
      geom_errorbar(aes(ymin=Mean_MAE-SD_MAE, 
                        ymax=Mean_MAE+SD_MAE), width=.2) +
      labs(title="Method Comparison: Mean Absolute Error",
           x="Method", y="MAE") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle=45, hjust=1))
    
    ggsave(file.path(output_dir, "method_comparison_mae.pdf"), p2, width=10, height=6)
  }
  
  # 3. Metrics by sample boxplot
  if (nrow(all_detailed_results) > 0) {
    p3 <- ggplot(all_detailed_results, aes(x=reorder(Method, PearsonCorr, FUN=median), y=PearsonCorr)) +
      geom_boxplot(fill="lightblue") +
      labs(title="Pearson Correlation by Method",
           x="Method", y="Pearson Correlation") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle=45, hjust=1))
    
    ggsave(file.path(output_dir, "sample_distribution_pearson.pdf"), p3, width=10, height=6)
  }
  
  cat("Results saved to:", output_dir, "\n")
  cat("Top performing methods by Pearson correlation:\n")
  print(head(ranked_methods[, c("Rank_Pearson", "Method", "Mean_PearsonCorr", "Mean_MAE", "SampleCount")], 5))
} else {
  cat("No valid results found to save\n")
}