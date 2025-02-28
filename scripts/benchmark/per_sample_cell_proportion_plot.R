#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
  stop(paste("5 arguments must be supplied instead of", length(args)), call. = FALSE)
}

####### Parameter of script (ORDER IS IMPORTANT)
path_Rlibrary <- args[1] #IMPORTANT
overall_ground_truth_path <- args[2]
per_sample_ground_truth_path <- args[3]
results_path <- args[4]
output_dir <- args[5]

# Libraries
.libPaths(path_Rlibrary, FALSE) #IMPORTANT
library(ggplot2)
library(reshape2)
library(dplyr)

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

# Load overall ground truth
load(overall_ground_truth_path)
overall_gt_proportions <- groundTruth$P

# Load per-sample ground truth
load(per_sample_ground_truth_path)
per_sample_gt_proportions <- groundTruth$P

# Get all unique cell types from both datasets (union instead of intersection)
all_cell_types <- unique(c(
  colnames(proportions),
  colnames(overall_gt_proportions),
  colnames(per_sample_gt_proportions)
))

# Print information about cell types for debugging
print(paste("Number of predicted cell types:", length(colnames(proportions))))
print(paste("Number of overall GT cell types:", length(colnames(overall_gt_proportions))))
print(paste("Number of per-sample GT cell types:", length(colnames(per_sample_gt_proportions))))
print(paste("Total unique cell types:", length(all_cell_types)))

# Check if there are any common cell types
common_cell_types <- Reduce(intersect, list(
  colnames(proportions),
  colnames(overall_gt_proportions),
  colnames(per_sample_gt_proportions)
))

if (length(common_cell_types) == 0) {
  print("WARNING: No common cell types found between ground truth and predicted proportions.")
  print("Will continue with the union of all cell types and fill missing values with zeros.")
}

# Helper function to expand a matrix with missing cell types (filled with zeros)
expand_matrix_with_zeros <- function(mat, all_cols) {
  current_cols <- colnames(mat)
  missing_cols <- setdiff(all_cols, current_cols)
  
  if (length(missing_cols) > 0) {
    # Create a matrix of zeros for the missing columns
    zero_mat <- matrix(0, nrow=nrow(mat), ncol=length(missing_cols))
    colnames(zero_mat) <- missing_cols
    rownames(zero_mat) <- rownames(mat)
    
    # Combine with the original matrix
    expanded_mat <- cbind(mat, zero_mat)
    # Reorder columns to match all_cols
    expanded_mat <- expanded_mat[, all_cols, drop=FALSE]
    return(expanded_mat)
  } else {
    # Reorder columns to match all_cols
    return(mat[, all_cols, drop=FALSE])
  }
}

# Expand all matrices to include all cell types
proportions <- expand_matrix_with_zeros(proportions, all_cell_types)
overall_gt_proportions <- expand_matrix_with_zeros(overall_gt_proportions, all_cell_types)
per_sample_gt_proportions <- expand_matrix_with_zeros(per_sample_gt_proportions, all_cell_types)

# Convert predicted proportions to long format
prop_df <- as.data.frame(proportions)
prop_df$Sample <- rownames(proportions)  # Ensure samples are correctly labeled
# Add a PairID column for grouping
prop_df$PairID <- prop_df$Sample
prop_long <- melt(prop_df, id.vars=c("Sample", "PairID"), 
                 variable.name="CellType", value.name="Proportion")
prop_long$Source <- "Predicted"
prop_long$DisplayOrder <- 1  # First in each pair

# Create overall ground truth data frame
overall_gt_df <- as.data.frame(overall_gt_proportions)
overall_gt_df$Sample <- "Overall GT"
overall_gt_df$PairID <- "ZZ_Overall" # Ensure it's at the end
overall_gt_long <- melt(overall_gt_df, id.vars=c("Sample", "PairID"), 
                       variable.name="CellType", value.name="Proportion")
overall_gt_long$Source <- "Ground Truth"
overall_gt_long$DisplayOrder <- 1  # Only one overall

# Create per-sample ground truth data frame
per_sample_gt_df <- as.data.frame(per_sample_gt_proportions)
per_sample_gt_df$Sample <- rownames(per_sample_gt_df)
# Match pair IDs with original samples
per_sample_gt_df$PairID <- per_sample_gt_df$Sample 
per_sample_gt_long <- melt(per_sample_gt_df, id.vars=c("Sample", "PairID"), 
                          variable.name="CellType", value.name="Proportion")
per_sample_gt_long$Source <- "Ground Truth"
per_sample_gt_long$DisplayOrder <- 2  # Second in each pair

# Get list of predicted samples and GT samples
predicted_samples <- unique(prop_long$Sample)
gt_samples <- unique(per_sample_gt_long$Sample)

# Find matching GT samples and non-matching GT samples
matching_gt_samples <- intersect(predicted_samples, gt_samples)
non_matching_gt_samples <- setdiff(gt_samples, predicted_samples)

# Add GT suffix to per-sample ground truth
per_sample_gt_long$Sample <- ifelse(
  per_sample_gt_long$Sample %in% matching_gt_samples,
  paste0(per_sample_gt_long$Sample, "-GT"),
  paste0("Unmatched-", per_sample_gt_long$Sample, "-GT")
)

# Exclude non-matching GT samples for simplicity
per_sample_gt_long <- per_sample_gt_long[!grepl("^Unmatched-", per_sample_gt_long$Sample), ]

# Combine all data (without final separator)
combined_data <- rbind(prop_long, per_sample_gt_long, overall_gt_long)

# Create a new column for actual display position with pair grouping
combined_data$DisplayPosition <- NA
current_pos <- 1

# For each unique pair ID (sorted)
for (pair in sort(unique(combined_data$PairID))) {
  # Skip overall which is handled separately
  if (pair %in% c("ZZ_Overall")) {
    next
  }
  
  # Get indices for this pair
  pair_indices <- which(combined_data$PairID == pair)
  
  # If there's only one sample in this pair (no matching GT)
  if (length(pair_indices) == length(all_cell_types)) {
    combined_data$DisplayPosition[pair_indices] <- current_pos
    current_pos <- current_pos + 1  # No extra spacing for unpaired samples
  } 
  # If there's a sample and its GT
  else {
    # Get display order within pair (1=sample, 2=GT)
    display_orders <- combined_data$DisplayOrder[pair_indices]
    
    # Assign positions within pair (close together)
    sample_indices <- pair_indices[display_orders == 1]
    gt_indices <- pair_indices[display_orders == 2]
    
    combined_data$DisplayPosition[sample_indices] <- current_pos
    combined_data$DisplayPosition[gt_indices] <- current_pos + 1
    
    current_pos <- current_pos + 3  # Skip one position for spacing between pairs
  }
}

# Handle overall GT position at the end
overall_indices <- which(combined_data$PairID == "ZZ_Overall")
combined_data$DisplayPosition[overall_indices] <- current_pos + 1  # Add extra spacing before Overall GT

# Add origin information to cell types
cell_types_in_gt <- unique(c(colnames(overall_gt_proportions), colnames(per_sample_gt_proportions)))
cell_types_in_pred <- colnames(proportions)

# Create a new column to mark the source of each cell type
combined_data$CellTypeSource <- "Both"
combined_data$CellTypeSource[combined_data$CellType %in% setdiff(cell_types_in_gt, cell_types_in_pred)] <- "GT only"
combined_data$CellTypeSource[combined_data$CellType %in% setdiff(cell_types_in_pred, cell_types_in_gt)] <- "Pred only"

# Add a suffix to cell type names to visually distinguish them
combined_data$CellTypeDisplay <- as.character(combined_data$CellType)
combined_data$CellTypeDisplay[combined_data$CellTypeSource == "GT only"] <- 
  paste0(combined_data$CellType[combined_data$CellTypeSource == "GT only"], " (GT)")
combined_data$CellTypeDisplay[combined_data$CellTypeSource == "Pred only"] <- 
  paste0(combined_data$CellType[combined_data$CellTypeSource == "Pred only"], " (Pred)")

# Define color palettes
gt_colors <- c(
  "#E57373",  # Light Red
  "#FF9933",  # Light Orange
  "#FFD54F",  # Light Gold
  "#AED581",  # Light Green
  "#4CBFFF",  # Light Cyan
  "#5999F2",  # Light Blue
  "#BA68C8",  # Light Purple
  "#E64CA6",  # Light Magenta
  "#FF8CBF",  # Light Pink
  "#7986CB"   # Light Navy Blue
)

# Darker colors for predicted-only cell types
pred_colors <- c(
  "#C62828",  # Dark Red
  "#E65100",  # Dark Orange
  "#F57F17",  # Dark Gold
  "#33691E",  # Dark Green
  "#01579B",  # Dark Cyan
  "#1A237E",  # Dark Blue
  "#4A148C",  # Dark Purple
  "#880E4F",  # Dark Magenta
  "#AD1457",  # Dark Pink
  "#283593"   # Dark Navy Blue
)
#### Wrong color assignment!!!!
# Get actual cell types from ground truth
gt_cell_types <- unique(c(colnames(overall_gt_proportions), colnames(per_sample_gt_proportions)))

# Create custom color mapping
unique_cell_displays <- unique(combined_data$CellTypeDisplay)
custom_colors <- character(length(unique_cell_displays))
names(custom_colors) <- unique_cell_displays

# Separate counters
light_counter <- 0
dark_counter <- 0

# Loop through each unique cell type display name
for (display_name in names(custom_colors)) {
  # Extract the base cell type name (remove any suffix like " (Pred)")
  base_cell_type <- gsub(" \\(Pred\\)| \\(GT\\)$", "", display_name)
  
  # If cell type appears in ground truth, use light colors
  if (base_cell_type %in% gt_cell_types) {
    light_counter <- light_counter + 1
    color_idx <- ((light_counter - 1) %% length(gt_colors)) + 1
    custom_colors[display_name] <- gt_colors[color_idx]
  } else {
    # Otherwise, it's prediction-only, use dark colors
    dark_counter <- dark_counter + 1
    color_idx <- ((dark_counter - 1) %% length(pred_colors)) + 1
    custom_colors[display_name] <- pred_colors[color_idx]
  }
}
###

# Print color assignments for debugging
print("Color assignments:")
for (i in 1:min(20, length(custom_colors))) {
  ct_display <- names(custom_colors)[i]
  row_idx <- which(cell_type_info$CellTypeDisplay == ct_display)[1]
  source <- cell_type_info$CellTypeSource[row_idx]
  print(paste(ct_display, "(", source, ") ->", custom_colors[i]))
}

# Create stacked bar plot
p <- ggplot(combined_data, aes(x=DisplayPosition, y=Proportion, fill=CellTypeDisplay)) +
  geom_bar(stat="identity", position="stack", width=0.8) +
  theme_minimal() +
  labs(title=paste(method_name, "vs Ground Truth Cell Type Proportions"),
       subtitle="Sample predictions paired with corresponding ground truth (non-matching cell types marked)",
       x="Samples", y="Proportions") +
  scale_x_continuous(
    breaks = combined_data$DisplayPosition,
    labels = combined_data$Sample,
    expand = c(0.02, 0)
  ) +
  theme(
        axis.text.x=element_text(angle=45, hjust=1),
        panel.grid.major=element_line(color="lightgray"),
        panel.grid.minor=element_line(color="lightgray"),
        legend.position="right",
        legend.key.size = unit(0.5, "cm"),
        plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust=0.5)
  ) +
  scale_fill_manual(values=custom_colors, name="Cell Type")

# Save the plot with adjusted dimensions based on the number of samples
n_samples <- length(unique(combined_data$Sample))
num_cell_types <- length(unique(combined_data$CellTypeDisplay))
plot_width <- max(12, n_samples * 0.4)  # Adjust width based on number of samples
plot_height <- 7 + (num_cell_types > 20) * (num_cell_types - 20) * 0.1  # Increase height for many cell types

output_filename <- file.path(output_dir, paste0(method_name, "_", data_name, "_paired_ground_truth.pdf"))
ggsave(output_filename, p, width=plot_width, height=plot_height)

print(paste("Plot saved to:", output_filename))