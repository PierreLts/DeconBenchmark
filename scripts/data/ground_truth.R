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

# Convert to long format for plotting - CORRECTED
# We need to transpose the matrix first to get samples as rows and cell types as columns
prop_df <- as.data.frame(t(proportions))
prop_df$Sample <- rownames(prop_df)
prop_long <- reshape2::melt(prop_df, id.vars="Sample", 
                            variable.name="CellType", value.name="Proportion")

# Create stacked bar plot
p <- ggplot(prop_long, aes(x=CellType, y=Proportion, fill=Sample)) +
  geom_bar(stat="identity", position="stack") +
  theme_minimal() +
  labs(title=paste(method_name, "Cell Type Proportions"),
       x="Samples", y="Proportions") +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        panel.grid.major=element_line(color="lightgray"),
        panel.grid.minor=element_line(color="lightgray"),
        legend.position="right",
        legend.key.size = unit(0.5, "cm"),
        plot.title=element_text(hjust=0.5)) +
  scale_fill_brewer(palette="Set3", name="Sample")

# Save the plot
output_filename <- file.path(output_dir, paste0(method_name, "_", data_name, "_cell_proportions.pdf"))
ggsave(output_filename, p, width=10, height=8)

print(paste("Plot saved to:", output_filename))