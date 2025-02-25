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

# Find all summary statistics files
summary_files <- list.files(
  path = benchmark_dir,
  pattern = "_summary_stats\\.csv$", 
  recursive = TRUE,
  full.names = TRUE
)

if (length(summary_files) == 0) {
  stop("No summary statistics files found in ", benchmark_dir)
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

# Reshape data for plotting
plot_data <- all_stats %>%
  filter(Metric %in% c("Pearson", "Spearman", "MAE", "RMSE", "R2", "CompositeScore")) %>%
  mutate(
    Label = sprintf("%.3f ± %.3f", Mean, SD),
    Lower = Mean - SD,
    Upper = Mean + SD
  )

# Rank models by Composite Score
method_order <- plot_data %>%
  filter(Metric == "CompositeScore") %>%
  arrange(desc(Mean)) %>%
  mutate(Rank = row_number()) %>%
  select(Method, Rank)

# Add ranking to method names
plot_data <- plot_data %>%
  left_join(method_order, by = "Method") %>%
  mutate(Method = paste0(Rank, ". ", Method))

# Update factor levels to maintain order
plot_data$Method <- factor(plot_data$Method, levels = unique(plot_data$Method[order(plot_data$Rank)]))

# Define metric display names
metric_display_names <- c(
  "Pearson (+)",
  "Spearman (+)",
  "MAE (-)",
  "RMSE (-)",
  "R² (+)",
  "Composite (+)"
)
names(metric_display_names) <- c("Pearson", "Spearman", "MAE", "RMSE", "R2", "CompositeScore")

# Grouped plot data
plot_data_grouped <- plot_data %>%
  mutate(
    MetricDisplay = factor(metric_display_names[Metric], levels = metric_display_names)
  )

# Find y-axis range
y_min <- min(plot_data_grouped$Lower, na.rm = TRUE)
y_max <- max(plot_data_grouped$Upper, na.rm = TRUE)
y_padding <- (y_max - y_min) * 0.1

# Define colors
metric_colors <- c(
  "Pearson (+)" = "#4CBFFF",
  "Spearman (+)" = "#AED581",
  "MAE (-)" = "#FF8CBF",
  "RMSE (-)" = "#FF9933",
  "R² (+)" = "#BA68C8",
  "Composite (+)" = "#E64CA6"
)

# Create single plot
single_plot <- ggplot(plot_data_grouped, aes(x = MetricDisplay, y = Mean, fill = MetricDisplay)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0, size = 0.7) + # Removed caps by setting width=0
  facet_grid(. ~ Method, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = metric_colors) +
  scale_y_continuous(limits = c(y_min - y_padding, y_max + y_padding)) +
  labs(
    title = "Performance Metrics Comparison Across Deconvolution Methods",
    subtitle = "(+) higher is better, (-) lower is better",
    y = "Score Value",
    x = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom",
    strip.text = element_text(size = 11, face = "bold"),
    panel.spacing.x = unit(0.25, "lines")
  ) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black")

# Save plot
single_file <- file.path(output_dir, "metrics_comparison.pdf")
ggsave(single_file, single_plot, width = max(12, length(unique(plot_data_grouped$Method)) * 1.5), height = 8)
print(paste("Saved metrics comparison plot to", single_file))

# Create summary table
summary_table <- all_stats %>%
  filter(Metric %in% names(metric_display_names)) %>%
  mutate(Value = sprintf("%.3f ± %.3f", Mean, SD)) %>%
  select(Method, Metric, Value) %>%
  pivot_wider(names_from = Metric, values_from = Value)

# Save summary table
summary_table_file <- file.path(output_dir, "all_methods_comparison.csv")
write_csv(summary_table, summary_table_file)
print(paste("Saved comparison table to", summary_table_file))
