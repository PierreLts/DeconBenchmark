#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop(paste("4 arguments must be supplied instead of", length(args)), call. = FALSE)
}

####### Parameter of script (ORDER IS IMPORTANT)
path_Rlibrary <- args[1] #IMPORTANT
ground_truth_path <- args[2]
results_path <- args[3]
output_dir <- args[4]

# Libraries
.libPaths(path_Rlibrary, FALSE) #IMPORTANT
library(ggplot2)
library(reshape2)

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

# Convert predicted proportions to long format (Samples on x-axis, Cell Types stacked)
prop_df <- as.data.frame(proportions)
prop_df$Sample <- rownames(proportions)  # Ensure samples are correctly labeled
prop_long <- melt(prop_df, id.vars="Sample", 
                 variable.name="CellType", value.name="Proportion")

# Load ground truth
load(ground_truth_path)
gt_proportions <- groundTruth$P

# Ensure ground truth has the same cell types as predicted values and reorder them
if (!all(colnames(proportions) %in% colnames(gt_proportions))) {
  stop("Mismatch: Ground truth and predicted proportions do not have the same cell types.")
}
gt_proportions <- gt_proportions[, colnames(proportions), drop=FALSE]  # Reorder columns to match

# Convert ground truth to a properly structured data frame
gt_df <- as.data.frame(gt_proportions)
gt_df$Sample <- "Ground Truth"  # Add identifier

# Introduce a small visual spacer by adding an empty label
spacer_df <- gt_df
spacer_df$Sample <- " "  # This will create a visual gap

# Convert ground truth to long format
gt_long <- melt(gt_df, id.vars="Sample", 
               variable.name="CellType", value.name="Proportion")

# Convert spacer to long format
spacer_long <- melt(spacer_df, id.vars="Sample", 
               variable.name="CellType", value.name="Proportion")
spacer_long$Proportion <- 0  # Ensure the spacer has no values

# Combine predicted, spacer, and ground truth data
combined_data <- rbind(prop_long, spacer_long, gt_long)

# Convert Sample to factor to maintain order (Spacer before Ground Truth)
combined_data$Sample <- factor(combined_data$Sample, levels = c(unique(prop_long$Sample), " ", "Ground Truth"))

# Generate a broader color palette for cell types
num_cell_types <- length(unique(combined_data$CellType))
custom_colors <- colorRampPalette(c(
  "#E57373",  # Light Red (lighter firebrick4)
  "#FF9933",  # Light Orange (lighter darkorange3)
  "#FFD54F",  # Light Gold (lighter gold3)
  "#AED581",  # Light Green (lighter chartreuse4)
  "#4CBFFF",  # Light Cyan (lighter deepskyblue3)
  "#5999F2",  # Light Blue (lighter steelblue4)
  "#BA68C8",  # Light Purple (lighter darkorchid4)
  "#E64CA6",  # Light Magenta (lighter mediumvioletred)
  "#FF8CBF",  # Light Pink (lighter hotpink3)
  "#7986CB"   # Light Navy Blue (lighter blue4)
)
)(num_cell_types)

# Assign colors to cell types
cell_types <- unique(combined_data$CellType)
names(custom_colors) <- cell_types  # Ensure colors are mapped correctly

# Create stacked bar plot (Samples on x-axis, Cell Types as fill)
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
  scale_fill_manual(values=custom_colors, name="Cell Type")  # Apply custom colors

# Save the plot
output_filename <- file.path(output_dir, paste0(method_name, "_", data_name, "_vs_ground_truth.pdf"))
ggsave(output_filename, p, width=12, height=6)

print(paste("Plot saved to:", output_filename))
