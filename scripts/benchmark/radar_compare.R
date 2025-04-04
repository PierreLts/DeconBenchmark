#!/usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript radar_compare.R <R_library_path> <output_file> <plot_title> <file1> <file2> ...")
}

# Set R library path
r_library_path <- args[1]
.libPaths(r_library_path, FALSE)

# Get output file path
original_output_file <- args[2]

# Get plot title
plot_title <- args[3]

# Get benchmark file paths
benchmark_files <- args[4:length(args)]

if (length(benchmark_files) == 0) {
  stop("No benchmark files provided")
}

# Get current working directory root (assuming we're in the project directory)
# This handles both absolute and relative paths based on current execution context
project_root <- getwd()
while (!dir.exists(file.path(project_root, "benchmark_results")) && 
       !file.exists(file.path(project_root, "DeconBenchmark")) && 
       project_root != dirname(project_root)) {
  project_root <- dirname(project_root)
}

# Construct target directory based on detected project structure
if (dir.exists(file.path(project_root, "DeconBenchmark", "benchmark_results"))) {
  target_dir <- file.path(project_root, "DeconBenchmark", "benchmark_results", "radar")
} else if (dir.exists(file.path(project_root, "benchmark_results"))) {
  target_dir <- file.path(project_root, "benchmark_results", "radar")
} else {
  # Fallback to relative path if we can't detect project structure
  target_dir <- file.path(getwd(), "benchmark_results", "radar")
}

cat("Using target directory:", target_dir, "\n")

# Create the directory if it doesn't exist
if (!dir.exists(target_dir)) {
  dir.create(target_dir, recursive = TRUE, showWarnings = TRUE)
  cat("Created output directory:", target_dir, "\n")
}

# Get the filename from the original path and create the new output path
output_filename <- basename(original_output_file)
output_file <- file.path(target_dir, output_filename)
cat("Output will be saved to:", output_file, "\n")

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
  
  # Replace "_benchmark_" with "_" to keep everything but remove "benchmark_"
  modified_name <- gsub("_benchmark_", "_", base_name)
  
  # Remove file extension
  modified_name <- gsub("\\.csv$", "", modified_name)
  
  return(modified_name)
}

# Function to parse HH:MM:SS time format to seconds
parse_time_to_seconds <- function(time_str) {
  # Ensure time_str is treated as character
  time_str <- as.character(time_str)
  
  # Debug output
  cat("  Parsing time string: '", time_str, "'\n", sep="")
  
  # Check if the string is already numeric
  if (!is.na(suppressWarnings(as.numeric(time_str)))) {
    cat("  Already numeric, returning as is:", as.numeric(time_str), "seconds\n")
    return(as.numeric(time_str))
  }
  
  # Handle N/A or empty strings
  if (is.na(time_str) || time_str == "" || tolower(time_str) == "n/a") {
    cat("  N/A or empty string detected\n")
    return(NA)
  }
  
  # Try to handle potential issues with time formatting
  time_str <- trimws(time_str)  # Remove any whitespace
  
  # Parse HH:MM:SS or MM:SS format
  parts <- strsplit(time_str, ":")[[1]]
  
  # Debug information about split parts
  cat("  Split into", length(parts), "parts:", paste(parts, collapse=", "), "\n")
  
  # Handle different time formats
  if (length(parts) == 3) {
    # HH:MM:SS format
    hours <- as.numeric(parts[1])
    minutes <- as.numeric(parts[2])
    seconds <- as.numeric(parts[3])
    total_seconds <- hours * 3600 + minutes * 60 + seconds
    cat("  Parsed as HH:MM:SS:", hours, "hours,", minutes, "minutes,", seconds, "seconds =", total_seconds, "total seconds\n")
    return(total_seconds)
  } else if (length(parts) == 2) {
    # MM:SS format
    minutes <- as.numeric(parts[1])
    seconds <- as.numeric(parts[2])
    total_seconds <- minutes * 60 + seconds
    cat("  Parsed as MM:SS:", minutes, "minutes,", seconds, "seconds =", total_seconds, "total seconds\n")
    return(total_seconds)
  } else if (length(parts) == 1) {
    # Try to interpret a single value - could be seconds or could be a different format
    if (grepl("^\\d+$", time_str)) {
      # If it's just a number, assume it's seconds
      seconds <- as.numeric(time_str)
      cat("  Parsed as seconds only:", seconds, "seconds\n")
      return(seconds)
    } else if (grepl("^\\d+s$", time_str, ignore.case = TRUE)) {
      # If it's in the format "123s", extract the number
      seconds <- as.numeric(sub("s$", "", time_str, ignore.case = TRUE))
      cat("  Parsed as seconds with 's' suffix:", seconds, "seconds\n")
      return(seconds)
    } else if (grepl("^\\d+m$", time_str, ignore.case = TRUE)) {
      # If it's in the format "123m", extract the number and convert to seconds
      minutes <- as.numeric(sub("m$", "", time_str, ignore.case = TRUE))
      seconds <- minutes * 60
      cat("  Parsed as minutes with 'm' suffix:", minutes, "minutes =", seconds, "seconds\n")
      return(seconds)
    }
  }
  
  # If we couldn't parse it using the above methods, try a more general approach
  # Check for patterns like "3m 20s" or "3 min 20 sec"
  if (grepl("\\d+\\s*m(in)?(ute)?s?", time_str, ignore.case = TRUE) && 
      grepl("\\d+\\s*s(ec)?(ond)?s?", time_str, ignore.case = TRUE)) {
    
    # Try to extract minutes and seconds from combined format
    minutes_match <- regexpr("\\d+\\s*m(in)?(ute)?s?", time_str, ignore.case = TRUE)
    minutes_str <- regmatches(time_str, minutes_match)
    minutes <- as.numeric(gsub("\\D", "", minutes_str))
    
    seconds_match <- regexpr("\\d+\\s*s(ec)?(ond)?s?", time_str, ignore.case = TRUE)
    seconds_str <- regmatches(time_str, seconds_match)
    seconds <- as.numeric(gsub("\\D", "", seconds_str))
    
    total_seconds <- minutes * 60 + seconds
    cat("  Parsed as combined format:", minutes, "minutes,", seconds, "seconds =", total_seconds, "seconds\n")
    return(total_seconds)
  }
  
  # Unable to parse
  warning("Unable to parse time format: ", time_str)
  cat("  Failed to parse time string\n")
  return(NA)
}

# Function to convert display values back to actual metric values
display_to_actual <- function(display_value, metric) {
  actual <- NA
  if (metric == "JSD") {
    actual <- 0.7 * (1 - display_value)
  } else if (metric == "R2") {
    actual <- (display_value * 11) - 10
  } else if (metric == "NRMSE") {
    actual <- 4 * (1 - display_value)
  } else if (metric == "PearsonCorr") {
    actual <- (display_value * 2) - 1
  } else if (metric == "Runtime") {
    actual <- 3600 * (1 - display_value)
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
    
    # Extract the overall values row (first row)
    metrics_row <- data[1, ]
    
    # Find columns for metrics
    # In the new format, we expect direct columns without "Models" prefix
    pearson_col <- which(colnames(data) %in% c("PearsonCorr", "Pearson"))[1]
    jsd_col <- which(colnames(data) == "JSD")[1]
    nrmse_col <- which(colnames(data) == "NRMSE")[1]
    r2_col <- which(colnames(data) %in% c("R2", "R^2"))[1]
    
    # For runtime, check exact column name first (prioritize "Runtime (H:M:S)")
    runtime_col <- which(colnames(data) == "Runtime (H:M:S)")[1]
    if (is.na(runtime_col)) {
      runtime_col <- which(colnames(data) %in% c("Runtime", "Time"))[1]
    }
    
    # Print what columns were found
    cat("  Found columns - Pearson:", 
        ifelse(!is.na(pearson_col), colnames(data)[pearson_col], "NOT FOUND"), 
        "JSD:", ifelse(!is.na(jsd_col), colnames(data)[jsd_col], "NOT FOUND"),
        "NRMSE:", ifelse(!is.na(nrmse_col), colnames(data)[nrmse_col], "NOT FOUND"),
        "R2:", ifelse(!is.na(r2_col), colnames(data)[r2_col], "NOT FOUND"),
        "Runtime:", ifelse(!is.na(runtime_col), colnames(data)[runtime_col], "NOT FOUND"), "\n")
    
    # If not found, try case-insensitive partial matches
    if (is.na(pearson_col)) pearson_col <- grep("pearson", colnames(data), ignore.case = TRUE)[1]
    if (is.na(jsd_col)) jsd_col <- grep("jsd", colnames(data), ignore.case = TRUE)[1]
    if (is.na(nrmse_col)) nrmse_col <- grep("nrmse", colnames(data), ignore.case = TRUE)[1]
    if (is.na(r2_col)) r2_col <- grep("r2|r\\^2", colnames(data), ignore.case = TRUE)[1]
    
    if (is.na(runtime_col)) {
      runtime_matches <- grep("runtime|time|h:m:s", colnames(data), ignore.case = TRUE)
      cat("  Potential runtime columns:", 
          ifelse(length(runtime_matches) > 0, 
                 paste(colnames(data)[runtime_matches], collapse=", "), 
                 "NONE"), "\n")
      runtime_col <- runtime_matches[1]  # Take the first match if multiple
    }
    
    # Check if all metrics were found
    if (any(is.na(c(pearson_col, jsd_col, nrmse_col, r2_col)))) {
      warning("Could not find all basic metrics in file: ", file_path)
      print(colnames(data))
      next
    }
    
    # Extract metrics values
    pearson_val <- as.numeric(metrics_row[[pearson_col]])
    jsd_val <- as.numeric(metrics_row[[jsd_col]])
    nrmse_val <- as.numeric(metrics_row[[nrmse_col]])
    r2_val <- as.numeric(metrics_row[[r2_col]])
    
    # Extract runtime value (if available)
    runtime_val <- NA
    if (!is.na(runtime_col)) {
      runtime_str <- metrics_row[[runtime_col]]
      cat("  Found runtime column:", colnames(data)[runtime_col], "with value:", runtime_str, "\n")
      runtime_val <- parse_time_to_seconds(runtime_str)
      cat("  Converted runtime to", runtime_val, "seconds\n")
    } else {
      warning("Runtime metric not found in file: ", file_path)
      cat("  Available columns:", paste(colnames(data), collapse=", "), "\n")
      runtime_val <- NA
    }
    
    cat("  Extracted metrics - Pearson:", pearson_val, "JSD:", jsd_val, "NRMSE:", nrmse_val, "R2:", r2_val, "Runtime:", runtime_val, "\n")
    
    # Put them in a data frame and add to our list
    data_list[[label]] <- data.frame(
      metric = c("PearsonCorr", "JSD", "NRMSE", "R2", "Runtime"),
      value = c(pearson_val, jsd_val, nrmse_val, r2_val, runtime_val),
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

# Check for minimum number of metrics
if (length(unique(all_data$metric)) < 3) {
  stop("At least 3 metrics are required for a meaningful radar plot")
}

# Check for single method case
if (length(unique(all_data$method)) == 1) {
  cat("Only one dataset provided - creating radar chart with single polygon\n")
}

# Transform values for consistent scale direction with specified ranges
all_data$display_value <- all_data$value  # Start with original values

# JSD: 0 (outer/best) to 0.7 (center/worst)
all_data$display_value[all_data$metric == "JSD"] <- 1 - pmin(1, all_data$value[all_data$metric == "JSD"] / 0.7)

# R2: 1 (outer/best) to -10 (center/worst)
all_data$display_value[all_data$metric == "R2"] <- (pmax(-10, all_data$value[all_data$metric == "R2"]) + 10) / 11

# NRMSE: 0 (outer/best) to 4 (center/worst)
all_data$display_value[all_data$metric == "NRMSE"] <- 1 - pmin(1, all_data$value[all_data$metric == "NRMSE"] / 4)

# Pearson Correlation: 1 (outer/best) to -1 (center/worst)
all_data$display_value[all_data$metric == "PearsonCorr"] <- (pmax(-1, pmin(1, all_data$value[all_data$metric == "PearsonCorr"])) + 1) / 2

# Runtime: 0 seconds (outer/best) to 3600 seconds (center/worst)
runtime_indices <- which(all_data$metric == "Runtime")
for (i in runtime_indices) {
  runtime_value <- all_data$value[i]
  method_name <- all_data$method[i]
  # Debug output for runtime values
  cat("Runtime value for", method_name, ":", runtime_value, "seconds\n")
  
  # Only transform if value is not NA
  if (!is.na(runtime_value)) {
    # Apply transformation, making sure 0 is best (1.0) and 3600+ is worst (0.0)
    display_value <- 1 - pmin(1, runtime_value / 3600)
    all_data$display_value[i] <- display_value
    cat("  Transformed to display value:", display_value, 
        "(1.0=fastest/0.0=slowest, based on 0-3600s scale)\n")
  } else {
    # For NA values, set to 0.5 (middle of the scale) to avoid missing data
    all_data$display_value[i] <- 0.5
    cat("  NA runtime value for", method_name, "- using default display value of 0.5\n")
  }
}

# Set factor level order for metrics to control display order
all_data$metric <- factor(all_data$metric, levels = c("PearsonCorr", "R2", "JSD", "NRMSE", "Runtime"))

# Set up colors - use the provided colors and continue the gradient
# Starting with your colors and continuing the gradient
colors <- c(
  "#FF7E1D",  # Orange
  "#DE0099",  # Pink
  "#B300F2",  # Purple
  "#6C13FF",  # Violet/Indigo
  "#0066FF",  # Blue
  "#00A3D7",  # Light Blue
  "#00C292",  # Teal
  "#00D858",  # Green
  "#78E100"   # Lime
)

# Create a proper color mapping regardless of the number of methods
methods <- unique(all_data$method)
num_methods <- length(methods)

if (num_methods > length(colors)) {
  # If more than 10 methods, create a colorful palette by interpolation
  warning("More than 10 methods provided. Using interpolated color palette.")
  colors <- colorRampPalette(colors)(num_methods)
}

# Assign colors to methods
color_mapping <- colors[1:num_methods]
names(color_mapping) <- methods

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
    label_data$y[i] <- label_data$y[i] + 0.1  # Reduced from 0.15 to 0.1
    label_data$x[i] <- label_data$x[i] - 0.03
  } else if (label_data$metric[i] == "Runtime") {
    # Adjust the Runtime label position if needed
    label_data$y[i] <- label_data$y[i] - 0.05
    label_data$x[i] <- label_data$x[i] - 0.05
  }
}

# Create data for scale value labels - only at key points (0, 0.5, 1)
scale_labels <- data.frame()

# Define universal offset for all metrics
default_offset <- 0.1  # Consistent shift for all metrics

for (i in seq_along(metrics)) {
  metric_name <- metrics[i]
  angle <- angles[i]
  
  # For each key break point (0, 0.5, 1)
  for (br in label_breaks) {
    # Calculate actual value for this display value
    actual <- display_to_actual(br, as.character(metric_name))
    
    # Calculate base position along the axis
    base_x <- br * cos(angle)
    base_y <- br * sin(angle)
    
    # Calculate position for the text label - aligned with axis direction
    # Use the angle to determine alignment
    hjust <- 0.5 # Default horizontal alignment
    vjust <- 0.5 # Default vertical alignment
    
    # Calculate offset to shift labels away from the axis
    # Offset direction depends on the axis angle
    if (angle == 0) { # Right
      hjust <- -0.1
      vjust <- 0.5
      x <- base_x + default_offset
      y <- base_y
    } else if (angle == pi/2) { # Top
      hjust <- 0.5
      vjust <- -0.1
      x <- base_x
      y <- base_y + default_offset
    } else if (angle == pi) { # Left
      hjust <- 1.1
      vjust <- 0.5
      x <- base_x - default_offset
      y <- base_y
    } else if (angle == 3*pi/2) { # Bottom
      hjust <- 0.5
      vjust <- 1.1
      x <- base_x
      y <- base_y - default_offset
    } else {
      # For diagonal axes, calculate the offset in both x and y directions
      offset_x <- default_offset * cos(angle)
      offset_y <- default_offset * sin(angle)
      x <- base_x + offset_x
      y <- base_y + offset_y
    }
    
    # Special adjustment for PearsonCorr
    if (metric_name == "PearsonCorr") {
      # Reduce the offset a bit for PearsonCorr to bring it closer
      adjustment_factor <- 0.7  # Reduce the offset to 70%
      x <- base_x + (x - base_x) * adjustment_factor
      y <- base_y + (y - base_y) * adjustment_factor
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
  scale_labels$metric == "Runtime",
  sprintf("%.0fs", scale_labels$actual_value),  # Format runtime in seconds
  ifelse(
    scale_labels$metric == "R2" | scale_labels$metric == "NRMSE", 
    sprintf("%.1f", scale_labels$actual_value),  # 1 decimal for R2 and NRMSE
    sprintf("%.2f", scale_labels$actual_value)   # 2 decimals for Pearson and JSD
  )
)

# Fixed plot size - no longer depends on number of methods
plot_size <- 10  # Fixed 10 inches regardless of dataset count

# Set up the plot
p <- ggplot() +
  # Add grid lines - CHANGED FROM grey85 TO black
  geom_path(
    data = grid_data %>% group_by(radius) %>% reframe(x = c(x, x[1]), y = c(y, y[1])),
    aes(x = x, y = y, group = radius),
    color = "black",
    size = 0.5
  ) +
  # Add radial axes - CHANGED FROM grey85 TO black
  geom_segment(
    data = data.frame(x = 0, y = 0, xend = axis_line_data$x, yend = axis_line_data$y),
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "black",
    size = 0.5
  ) +
  # Add lines connecting points (using path with Cartesian coordinates)
  geom_polygon(
    data = polygon_data %>% group_by(method) %>% reframe(
      x = c(x, x[1]),
      y = c(y, y[1])
    ),
    aes(x = x, y = y, group = method, color = method),
    fill = NA,
    size = 0.8
  ) +
  # Add the points
  geom_point(
    data = polygon_data,
    aes(x = x, y = y, color = method),
    size = 2
  ) +
  # Add scale value labels - CHANGED FROM grey40 TO black
  geom_text(
    data = scale_labels,
    aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
    color = "black",
    size = 3.5
  ) +
  # Use Cartesian coordinates with equal aspect ratio
  coord_equal(xlim = c(-2, 2), ylim = c(-1.6, 1.6), expand = FALSE) +
  # Use custom colors
  scale_color_manual(values = color_mapping) +
  # Add title
  labs(title = plot_title) +
  # Add axis labels - use label_data for adjusted positions
  geom_text(
    data = label_data,
    aes(x = x, y = y, label = metric),
    hjust = 0.5,
    vjust = 0.5,
    size = 5,
    fontface = "bold"
  )

# Fixed: Always use 10 datasets per column for legend
legend_items_per_column <- 10  # Always 10 datasets per column 
legend_columns <- ceiling(num_methods / legend_items_per_column)

# Position the legend in the top right corner
legend_x <- 0.8  # Position on the right side

# Reduce spacing between legend items for many datasets
legend_spacing <- 0.08  # Tighter spacing for many datasets

col_width <- 0.5  # Width of each column
legend_start_y <- 1.3  # Starting position

# Create the legend data
legend_data <- data.frame(
  x = numeric(),
  y = numeric(),
  label = character(),
  method = character(),
  stringsAsFactors = FALSE
)

# Calculate positions for each legend item
for (col in 1:legend_columns) {
  x_pos <- legend_x + (col - 1) * col_width
  
  # Calculate starting and ending indices for this column
  start_idx <- (col - 1) * legend_items_per_column + 1
  end_idx <- min(col * legend_items_per_column, length(methods))
  
  # Calculate y positions for this column's items
  y_positions <- seq(from = legend_start_y, by = -legend_spacing, length.out = end_idx - start_idx + 1)
  
  for (i in start_idx:end_idx) {
    item_idx <- i - start_idx + 1
    method_name <- methods[i]
    
    legend_data <- rbind(legend_data, data.frame(
      x = x_pos,
      y = y_positions[item_idx],
      label = method_name,
      method = method_name,
      stringsAsFactors = FALSE
    ))
  }
}

# Create theme settings
theme_settings <- theme_void() +
  theme(
    # Center title
    plot.title = element_text(hjust = 0.5, size = 16),
    # Hide the default legend
    legend.position = "none",
    # Background and margins
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(20, 20, 20, 20)
  )

# Apply theme settings
p <- p + theme_settings

# Adjust legend font size based on number of methods
legend_font_size <- 4.0 

# Add colored text annotations as legend
p <- p + 
  geom_text(
    data = legend_data,
    aes(x = x, y = y, label = label, color = method),
    size = legend_font_size,
    fontface = "bold",
    hjust = 0,  # Left-align the text
    vjust = 0.5  # Center vertically
  )

# Scale down width adjustment for many datasets to prevent excessive width
legend_width_adjustment <- 0.8

# Cap the maximum width to ensure it doesn't get too wide
max_width <- 15  # Maximum width in inches
final_width <- min(plot_size + legend_width_adjustment, max_width)

# Save plot with controlled dimensions - add create.dir=TRUE
ggsave(output_file, p, width = final_width, height = plot_size, units = "in", create.dir = TRUE)
cat("Radar plot saved to:", output_file, "\n")