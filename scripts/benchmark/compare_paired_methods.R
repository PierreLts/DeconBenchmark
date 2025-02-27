#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
  stop("3 arguments must be supplied: <R_library_path> <benchmark_results_dir> <output_dir>", call.=FALSE)
}

####### Parameters
path_Rlibrary <- args[1]
benchmark_dir <- args[2]
output_dir <- args[3]

# Libraries
.libPaths(path_Rlibrary, FALSE)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Find all paired summary statistics files
summary_files <- list.files(
  path = benchmark_dir,
  pattern = "_paired_summary_stats\\.csv$", 
  recursive = TRUE,
  full.names = TRUE
)

if (length(summary_files) == 0) {
  stop("No paired summary statistics files found in ", benchmark_dir)
}

# Function to extract method name from file path
extract_method_name <- function(file_path) {
  filename <- basename(file_path)
  method_name <- str_extract(filename, "^[^_]+")
  return(method_name)
}

# Read and combine all summary statistics
all_stats <- data.frame()
for (file in summary_files) {
  method <- extract_method_name(file)
  stats <- read_csv(file)
  
  stats$Method <- method
  
  all_stats <- bind_rows(all_stats, stats)
}

# Ensure numeric columns are numeric
numeric_cols <- c("Mean", "SD")
all_stats[numeric_cols] <- lapply(all_stats[numeric_cols], as.numeric)

# Filter for metrics of interest
metrics_of_interest <- c("Pearson", "Spearman", "MAE", "RMSE", "R2", "CompositeScore")
plot_data <- all_stats %>%
  filter(Metric %in% metrics_of_interest)

# Create performance comparison for each type of comparison (Paired, Unpaired, Overall)
for (comparison_type in unique(plot_data$Type)) {
  # Filter for this comparison type
  type_data <- plot_data %>%
    filter(Type == comparison_type)
  
  # Create a separate summary for each metric
  for (metric_name in metrics_of_interest) {
    metric_data <- type_data %>%
      filter(Metric == metric_name)
    
    # Determine if higher is better for this metric
    higher_better <- !(metric_name %in% c("MAE", "RMSE"))
    
    # Order methods by this metric
    if (higher_better) {
      ordered_methods <- metric_data %>%
        arrange(desc(Mean)) %>%
        pull(Method)
    } else {
      ordered_methods <- metric_data %>%
        arrange(Mean) %>%
        pull(Method)
    }
    
    # Add ranking
    metric_data <- metric_data %>%
      mutate(
        Method = factor(Method, levels = ordered_methods),
        Rank = row_number(),
        Label = sprintf("%.3f ± %.3f", Mean, SD),
        Direction = if(higher_better) "(+)" else "(-)"
      )
    
    # Set plot title
    plot_title <- paste0(
      metric_name, " ", 
      if(higher_better) "(higher is better)" else "(lower is better)",
      " - ", comparison_type, " Comparison"
    )
    
    # Create the plot
    p <- ggplot(metric_data, aes(x = reorder(Method, if(higher_better) Mean else -Mean), y = Mean)) +
      geom_bar(stat = "identity", fill = if(higher_better) "#4CBFFF" else "#FF8CBF") +
      geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2) +
      labs(
        title = plot_title,
        subtitle = paste0("Methods ranked by ", metric_name, " value"),
        x = "Method",
        y = paste0(metric_name, " Value")
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)
      ) +
      geom_text(aes(label = Label, y = Mean + SD + if(higher_better) 0.05 else -0.05),
                hjust = if(higher_better) 0 else 1,
                angle = 90, 
                size = 3)
    
    # Save the plot
    metric_file_name <- gsub("[^a-zA-Z0-9]", "_", metric_name)
    type_file_name <- gsub("[^a-zA-Z0-9]", "_", comparison_type)
    output_filename <- file.path(output_dir, paste0(metric_file_name, "_", type_file_name, "_comparison.pdf"))
    ggsave(output_filename, p, width = max(10, length(unique(metric_data$Method)) * 0.5), height = 8)
    print(paste("Saved", metric_name, "comparison for", comparison_type, "to:", output_filename))
  }
  
  # Create combined visualization for this comparison type
  type_data_wide <- type_data %>%
    pivot_wider(
      id_cols = c(Method),
      names_from = Metric,
      values_from = Mean
    )
  
  # Rank methods by composite score
  type_data_wide <- type_data_wide %>%
    arrange(desc(CompositeScore)) %>%
    mutate(
      Method_Ranked = paste0(row_number(), ". ", Method),
      Method_Ranked = factor(Method_Ranked, levels = Method_Ranked)
    )
  
  # Calculate rankings for each metric
  ranked_data <- type_data %>%
    group_by(Metric) %>%
    mutate(
      Rank = if_else(
        Metric %in% c("MAE", "RMSE"),
        dense_rank(Mean),
        dense_rank(desc(Mean))
      )
    ) %>%
    ungroup() %>%
    select(Method, Metric, Rank) %>%
    pivot_wider(
      id_cols = Method,
      names_from = Metric,
      values_from = Rank,
      names_prefix = "Rank_"
    ) %>%
    right_join(type_data_wide, by = "Method")
  
  # Create a heatmap of scores
  plot_data_long <- ranked_data %>%
    pivot_longer(
      cols = starts_with("Rank_"),
      names_to = "Metric",
      values_to = "Rank"
    ) %>%
    mutate(
      Metric = gsub("Rank_", "", Metric),
      Metric = factor(Metric, levels = metrics_of_interest)
    )
  
  # Create heatmap
  heatmap_plot <- ggplot(plot_data_long, aes(x = Metric, y = Method_Ranked, fill = Rank)) +
    geom_tile() +
    scale_fill_gradient(low = "#4CBFFF", high = "#FF8CBF", name = "Rank") +
    labs(
      title = paste("Method Rankings by Metric -", comparison_type, "Comparison"),
      subtitle = "Lower rank (blue) is better",
      x = "Metric",
      y = "Method"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  # Save the heatmap
  type_file_name <- gsub("[^a-zA-Z0-9]", "_", comparison_type)
  heatmap_filename <- file.path(output_dir, paste0("ranking_heatmap_", type_file_name, ".pdf"))
  ggsave(heatmap_filename, heatmap_plot, width = 10, height = max(8, nrow(ranked_data) * 0.4))
  print(paste("Saved ranking heatmap for", comparison_type, "to:", heatmap_filename))
  
  # Create a table with all metrics
  table_data <- ranked_data %>%
    select(-starts_with("Rank_")) %>%
    mutate_if(is.numeric, round, 3)
  
  # Save as CSV
  table_filename <- file.path(output_dir, paste0("all_metrics_", type_file_name, ".csv"))
  write_csv(table_data, table_filename)
  print(paste("Saved metrics table for", comparison_type, "to:", table_filename))
}

# Create a combined table with all comparison types
all_types_summary <- all_stats %>%
  mutate(
    Value = sprintf("%.3f ± %.3f", Mean, SD)
  ) %>%
  select(Method, Metric, Type, Value) %>%
  pivot_wider(
    id_cols = c(Method, Metric),
    names_from = Type,
    values_from = Value
  )

# Save the combined table
combined_table_filename <- file.path(output_dir, "all_methods_all_types_comparison.csv")
write_csv(all_types_summary, combined_table_filename)
print(paste("Saved combined comparison table to:", combined_table_filename))

print("Comparison of paired benchmarking results complete.")