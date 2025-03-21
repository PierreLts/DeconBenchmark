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

# Function to convert display values back to actual metric values
display_to_actual <- function(display_value, metric) {
  actual <- NA
  if (metric == "JSD") {
    actual <- 0.5 * (1 - display_value)
  } else if (metric == "R2") {
    actual <- (display_value * 16) - 15
  } else if (metric == "NRMSE") {
    actual <- 5 * (1 - display_value)
  } else if (metric == "PearsonCorr") {
    actual <- (display_value * 2) - 1
  }
  return(actual)
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
    if (nrow(data) < 1) {
      warning("File has fewer than 2 rows: ", file_path)
      next
    }
    
    # Extract the overall values row (second row)
    metrics_row <- data[1, ]
    
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
    r2_val <- as.numeric(metrics_row[[r2_col]])
    
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

# Transform values for consistent scale direction with specified ranges
all_data$display_value <- all_data$value  # Start with original values

# JSD: 0 (outer/best) to 0.2 (center/worst)
all_data$display_value[all_data$metric == "JSD"] <- 1 - pmin(1, all_data$value[all_data$metric == "JSD"] / 0.5)

# R2: 1 (outer/best) to -4 (center/worst)
all_data$display_value[all_data$metric == "R2"] <- (pmax(-15, all_data$value[all_data$metric == "R2"]) + 15) / 16

# NRMSE: 0 (outer/best) to 5 (center/worst)
all_data$display_value[all_data$metric == "NRMSE"] <- 1 - pmin(1, all_data$value[all_data$metric == "NRMSE"] / 5)

# Pearson Correlation: 1 (outer/best) to -1 (center/worst)
all_data$display_value[all_data$metric == "PearsonCorr"] <- (pmax(-1, pmin(1, all_data$value[all_data$metric == "PearsonCorr"])) + 1) / 2

# Set factor level order for metrics to control display order
all_data$metric <- factor(all_data$metric, levels = c("PearsonCorr", "R2", "JSD", "NRMSE"))

# Set up colors - use distinctive colors
colors <- c(
  "#1E90FF",  # Dodger Blue - for select-A
  "#32CD32",  # Lime Green - for select-AB
  "#FF69B4"   # Hot Pink - for select-B
)
names(colors) <- unique(all_data$method)

# Function to calculate Cartesian coordinates for true straight lines in radar plot
create_polygon_data <- function(df) {
  methods <- unique(df$method)
  metrics <- levels(df$metric)
  n_metrics <- length(metrics)
  
  result <- data.frame()
  
  for (m in methods) {
    # Filter data for this method
    method_data <- df[df$method == m, ]
    method_data <- method_data[order(method_data$metric), ]
    
    # Convert to polar coordinates for plotting
    angles <- (as.numeric(method_data$metric) - 1) * 2 * pi / n_metrics
    radii <- method_data$display_value
    
    # Convert to Cartesian coordinates
    x <- radii * cos(angles)
    y <- radii * sin(angles)
    
    # Create a data frame for this method
    method_df <- data.frame(
      x = x,
      y = y,
      radius = radii,
      angle = angles,
      method = m,
      metric = method_data$metric
    )
    
    result <- rbind(result, method_df)
  }
  
  return(result)
}

# Create data for plotting
polygon_data <- create_polygon_data(all_data)

# Define the breaks for the radial grid lines - keep all intermediate lines
grid_breaks <- c(0, 0.25, 0.5, 0.75, 1)

# But only show labels at key points
label_breaks <- c(0, 0.5, 1)

# Prepare grid lines data for radar chart
metrics <- levels(all_data$metric)
n_metrics <- length(metrics)
angles <- (seq_along(metrics) - 1) * 2 * pi / n_metrics

grid_data <- expand.grid(
  angle = angles,
  radius = grid_breaks
)

grid_data$x <- grid_data$radius * cos(grid_data$angle)
grid_data$y <- grid_data$radius * sin(grid_data$angle)

# Create axis line data (for the straight lines)
axis_line_data <- data.frame(
  angle = angles,
  radius = 1.0,  # Extend to edge of radar plot
  x = 1.0 * cos(angles),
  y = 1.0 * sin(angles),
  metric = metrics
)

# Create separate data for labels that can be adjusted independently
label_data <- data.frame(
  angle = angles,
  radius = 1.25,  # Base distance from chart edge
  x = 1.25 * cos(angles),
  y = 1.25 * sin(angles),
  metric = metrics
)

# Adjust only the label positions, not the axis lines
for (i in 1:nrow(label_data)) {
  if (label_data$metric[i] == "JSD") {
    label_data$y[i] <- label_data$y[i] + 0.15
    label_data$x[i] <- label_data$x[i] + 0.03
  } else if (label_data$metric[i] == "PearsonCorr") {
    label_data$y[i] <- label_data$y[i] + 0.15
    label_data$x[i] <- label_data$x[i] - 0.03
  }
}

# Create data for scale value labels - only at key points (0, 0.5, 1)
scale_labels <- data.frame()

for (i in seq_along(metrics)) {
  metric_name <- metrics[i]
  angle <- angles[i]
  
  # For each key break point (0, 0.5, 1)
  for (br in label_breaks) {
    # Calculate actual value for this display value
    actual <- display_to_actual(br, as.character(metric_name))
    
    # Calculate label position along the axis
    x <- br * cos(angle)
    y <- br * sin(angle)
    
    # Calculate position for the text label - aligned with axis direction
    # Use the angle to determine alignment
    hjust <- 0.5 # Default horizontal alignment
    vjust <- 0.5 # Default vertical alignment
    
    # Special case for center (0,0)
    if (br == 0) {
      # For center value, slightly offset in the direction of the axis
      # to avoid overlapping at the center
      center_offset <- 0.03  # Reduced from 0.05 to be closer to center
      x <- center_offset * cos(angle)
      y <- center_offset * sin(angle)
    }
    
    # Adjust alignment based on axis quadrant
    if (angle == 0) { # Right
      hjust <- -0.1
      vjust <- 0.5
      x <- x + 0.05
    } else if (angle == pi/2) { # Top
      hjust <- 0.5
      vjust <- -0.1
      y <- y + 0.05
    } else if (angle == pi) { # Left
      hjust <- 1.1
      vjust <- 0.5
      x <- x - 0.05
    } else if (angle == 3*pi/2) { # Bottom
      hjust <- 0.5
      vjust <- 1.1
      y <- y - 0.05
    }
    
    # Add to data frame
    scale_labels <- rbind(scale_labels, data.frame(
      metric = metric_name,
      display_value = br,
      actual_value = actual,
      x = x,
      y = y,
      hjust = hjust,
      vjust = vjust
    ))
  }
}

# Format scale labels (round to appropriate decimal places)
scale_labels$label <- ifelse(
  scale_labels$metric == "R2" | scale_labels$metric == "NRMSE", 
  sprintf("%.1f", scale_labels$actual_value),  # 1 decimal for R2 and NRMSE
  sprintf("%.2f", scale_labels$actual_value)   # 2 decimals for Pearson and JSD
)

# Set up the plot
p <- ggplot() +
  # Add grid lines
  geom_path(
    data = grid_data %>% group_by(radius) %>% summarize(x = c(x, x[1]), y = c(y, y[1])),
    aes(x = x, y = y, group = radius),
    color = "grey85",
    size = 0.5
  ) +
  # Add radial axes - use axis_line_data for straight lines
  geom_segment(
    data = data.frame(x = 0, y = 0, xend = axis_line_data$x, yend = axis_line_data$y),
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "grey85",
    size = 0.5
  ) +
  # Add lines connecting points (using path with Cartesian coordinates)
  geom_polygon(
    data = polygon_data %>% group_by(method) %>% summarize(
      x = c(x, x[1]),
      y = c(y, y[1])
    ),
    aes(x = x, y = y, group = method, color = method),
    fill = NA,
    size = 1.2
  ) +
  # Add the points
  geom_point(
    data = polygon_data,
    aes(x = x, y = y, color = method),
    size = 3
  ) +
  # Add scale value labels with proper alignment
  geom_text(
    data = scale_labels,
    aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
    color = "grey40",
    size = 3
  ) +
  # Use Cartesian coordinates with equal aspect ratio
  coord_equal() +
  # Use custom colors
  scale_color_manual(values = colors) +
  # Add title
  labs(title = "Performance Metrics Comparison") +
  # Add axis labels - use label_data for adjusted positions
  geom_text(
    data = label_data,
    aes(x = x, y = y, label = metric),
    hjust = 0.5,
    vjust = 0.5,
    size = 4,
    fontface = "bold"
  ) +
  # Set theme
  theme_void() +
  theme(
    # Center title
    plot.title = element_text(hjust = 0.5, size = 16),
    # Format legend
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    # Background and margins
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(20, 20, 20, 20)
  )

# Save plot
ggsave(output_file, p, width = 8, height = 8, units = "in")
cat("Radar plot saved to:", output_file, "\n")