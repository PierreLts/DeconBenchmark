#!/usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript radar_compare.R <R_library_path> <output_file> <file1> <file2> ...")
}

# Set R library path
r_library_path <- args[1]
.libPaths(r_library_path, FALSE)

# Get output file path
output_file <- args[2]

# Get benchmark file paths
benchmark_files <- args[3:length(args)]

# Function to install packages if not available
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing package:", pkg, "\n")
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

# Install and load required packages
install_if_missing("ggplot2")
install_if_missing("readr")
install_if_missing("dplyr")
install_if_missing("tidyr")
install_if_missing("scales")

# Function to extract dataset name and selection from file path
extract_name_and_selection <- function(file_path) {
  base_name <- basename(file_path)
  
  # Extract dataset name (everything before "_benchmark_")
  dataset_name <- gsub("_benchmark_.*$", "", base_name)
  
  # Extract selection (if present)
  if (grepl("select-", base_name)) {
    selection <- gsub(".*_select-([^\\.]+)\\.csv$", "\\1", base_name)
    return(paste0(dataset_name, "_select-", selection))
  }
  
  return(dataset_name)
}

# Read and process benchmark files
data_list <- list()
labels <- character()

for (file_path in benchmark_files) {
  # Extract dataset name and selection for labeling
  label <- extract_name_and_selection(file_path)
  labels <- c(labels, label)
  
  cat("Processing file:", file_path, "\n")
  
  # Read CSV file
  tryCatch({
    data <- read_csv(file_path)
    
    # Check file structure
    if (nrow(data) < 2) {
      warning("File has fewer than 2 rows: ", file_path)
      next
    }
    
    # Extract the overall values row (second row)
    metrics_row <- data[2, ]
    
    # Find columns for metrics
    # First try exact matches
    pearson_col <- which(colnames(data) %in% c("PearsonCorr", "Pearson"))[1]
    jsd_col <- which(colnames(data) == "JSD")[1]
    nrmse_col <- which(colnames(data) == "NRMSE")[1]
    r2_col <- which(colnames(data) %in% c("R2", "R^2"))[1]
    
    # If not found, try case-insensitive partial matches
    if (is.na(pearson_col)) pearson_col <- grep("pearson", colnames(data), ignore.case = TRUE)[1]
    if (is.na(jsd_col)) jsd_col <- grep("jsd", colnames(data), ignore.case = TRUE)[1]
    if (is.na(nrmse_col)) nrmse_col <- grep("nrmse", colnames(data), ignore.case = TRUE)[1]
    if (is.na(r2_col)) r2_col <- grep("r2|r\\^2", colnames(data), ignore.case = TRUE)[1]
    
    # Check if all metrics were found
    if (any(is.na(c(pearson_col, jsd_col, nrmse_col, r2_col)))) {
      warning("Could not find all metrics in file: ", file_path)
      print(colnames(data))
      next
    }
    
    # Extract metrics values
    pearson_val <- as.numeric(metrics_row[[pearson_col]])
    jsd_val <- as.numeric(metrics_row[[jsd_col]])
    nrmse_val <- as.numeric(metrics_row[[nrmse_col]])
    r2_val <- as.numeric(metrics_row[[r2_col]])  # FIXED: Was using r2_val instead of r2_col
    
    cat("  Extracted metrics - Pearson:", pearson_val, "JSD:", jsd_val, "NRMSE:", nrmse_val, "R2:", r2_val, "\n")
    
    # Put them in a data frame and add to our list
    data_list[[label]] <- data.frame(
      metric = c("PearsonCorr", "JSD", "NRMSE", "R2"),
      value = c(pearson_val, jsd_val, nrmse_val, r2_val),
      method = label
    )
    
  }, error = function(e) {
    warning("Error processing file ", file_path, ": ", e$message)
  })
}

if (length(data_list) == 0) {
  stop("No valid data found in any of the provided files")
}

# Combine all data frames into one
all_data <- do.call(rbind, data_list)

# Transform values for consistent scale direction
# For JSD and NRMSE, lower is better, so we invert them
all_data$display_value <- all_data$value
all_data$display_value[all_data$metric == "JSD"] <- 1 - all_data$value[all_data$metric == "JSD"]
all_data$display_value[all_data$metric == "NRMSE"] <- 1 - (all_data$value[all_data$metric == "NRMSE"] / 
                                                          max(all_data$value[all_data$metric == "NRMSE"]))

# Set factor level order for metrics to control display order
all_data$metric <- factor(all_data$metric, levels = c("PearsonCorr", "R2", "JSD", "NRMSE"))

# Instead of using geom_path with polar coordinates which creates curved lines,
# we'll manually create straight line segments between points in polar coordinates
create_straight_segments <- function(df) {
  metrics <- levels(df$metric)
  methods <- unique(df$method)
  n_metrics <- length(metrics)
  
  # Create a list to store segment data frames
  segment_dfs <- list()
  
  for (m in methods) {
    method_data <- df[df$method == m, ]
    
    # Make sure data is ordered by metric factor levels
    method_data <- method_data[order(method_data$metric), ]
    
    # Add the first point at the end to close the polygon
    method_data <- rbind(method_data, method_data[1, ])
    
    # Create segment data with x1,y1,x2,y2 coordinates for each pair of points
    for (i in 1:(nrow(method_data) - 1)) {
      # Get adjacent points
      p1 <- method_data[i, ]
      p2 <- method_data[i + 1, ]
      
      # Convert metrics to angles (in radians)
      x1 <- (match(as.character(p1$metric), metrics) - 1) * 2 * pi / n_metrics
      x2 <- (match(as.character(p2$metric), metrics) - 1) * 2 * pi / n_metrics
      
      # For the closing segment, adjust angle
      if (i == n_metrics) {
        x2 <- 2 * pi  # Return to the starting point
      }
      
      # Use display_value for radius
      y1 <- p1$display_value
      y2 <- p2$display_value
      
      # Store as a segment
      segment_dfs[[length(segment_dfs) + 1]] <- data.frame(
        x1 = x1, y1 = y1, 
        x2 = x2, y2 = y2,
        method = m
      )
    }
  }
  
  # Combine all segments
  do.call(rbind, segment_dfs)
}

# Create data with straight line segments
segments_data <- create_straight_segments(all_data)




# Set up colors - use distinctive colors
colors <- c(
  "#1E90FF",  # Dodger Blue - for select-A
  "#32CD32",  # Lime Green - for select-AB
  "#FF69B4"   # Hot Pink - for select-B
)
names(colors) <- unique(all_data$method)

# Instead of using geom_path with polar coordinates which creates curved lines,
# we'll manually create straight line segments between points in polar coordinates
create_straight_segments <- function(df) {
  metrics <- levels(df$metric)
  methods <- unique(df$method)
  n_metrics <- length(metrics)
  
  # Create a list to store segment data frames
  segment_dfs <- list()
  
  for (m in methods) {
    method_data <- df[df$method == m, ]
    
    # Make sure data is ordered by metric factor levels
    method_data <- method_data[order(method_data$metric), ]
    
    # Add the first point at the end to close the polygon
    method_data <- rbind(method_data, method_data[1, ])
    
    # Create segment data with x1,y1,x2,y2 coordinates for each pair of points
    for (i in 1:(nrow(method_data) - 1)) {
      # Get adjacent points
      p1 <- method_data[i, ]
      p2 <- method_data[i + 1, ]
      
      # Convert metrics to angles (in radians)
      x1 <- (match(as.character(p1$metric), metrics) - 1) * 2 * pi / n_metrics
      x2 <- (match(as.character(p2$metric), metrics) - 1) * 2 * pi / n_metrics
      
      # For the closing segment, adjust angle
      if (i == n_metrics) {
        x2 <- 2 * pi  # Return to the starting point
      }
      
      # Use display_value for radius
      y1 <- p1$display_value
      y2 <- p2$display_value
      
      # Store as a segment
      segment_dfs[[length(segment_dfs) + 1]] <- data.frame(
        x1 = x1, y1 = y1, 
        x2 = x2, y2 = y2,
        method = m
      )
    }
  }
  
  # Combine all segments
  do.call(rbind, segment_dfs)
}

# Create data with straight line segments
segments_data <- create_straight_segments(all_data)

# Define the breaks for the radial grid lines
breaks <- c(0, 0.25, 0.5, 0.75, 1)

# Set up the plot
# Set up the plot
p <- ggplot() +
  # Add the straight line segments
  geom_segment(
    data = segments_data,
    aes(
      x = x1, y = y1, 
      xend = x2, yend = y2,
      color = method
    ),
    linewidth = 1.2
  ) +
  # Add the points on top
  geom_point(
    data = all_data,
    aes(
      x = (as.numeric(metric) - 1) * 2 * pi / length(levels(metric)),
      y = display_value,
      color = method
    ),
    size = 3
  ) +
  # Use polar coordinates with theta starting at top (0=north, pi/2=east)
  coord_polar(start = pi/2) +
  # Use the same colors for lines and points
  scale_color_manual(values = colors) +
  # Add title
  labs(title = "Performance Metrics Comparison") +
  # Set y-axis breaks for the grid circles
  scale_y_continuous(
    breaks = breaks,
    labels = breaks,
    limits = c(0, 1),
    expand = expansion(mult = c(0.1, 0.2))  # Expand the plot margins
  ) +
  # Add axis labels at the correct positions
  scale_x_continuous(
    breaks = seq(0, 2*pi, length.out = length(levels(all_data$metric)) + 1)[1:4],
    labels = levels(all_data$metric),
    limits = c(0, 2*pi)
  ) +
  # Remove default axis labels
  theme_minimal() +
  theme(
    # Center title
    plot.title = element_text(hjust = 0.5, size = 16),
    # Format axis text (metric labels)
    axis.text.x = element_text(size = 12, margin = margin(t = 15, unit = "pt")),
    axis.text.y = element_text(size = 9),
    axis.title = element_blank(),
    # Format legend
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    # Format grid
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_line(color = "grey95")
  )

# Save plot
ggsave(output_file, p, width = 8, height = 8, units = "in")
cat("Radar plot saved to:", output_file, "\n")