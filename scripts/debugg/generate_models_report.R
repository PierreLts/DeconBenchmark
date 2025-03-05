#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
  stop("2 arguments must be supplied: R_LIBRARY_PATH OUTPUT_DIRECTORY", call.=FALSE)
}

# Set library path
path_Rlibrary <- args[1]
output_dir <- args[2]
.libPaths(path_Rlibrary, FALSE)

# Load library
library(DeconBenchmark)
output_file <- file.path(output_dir, "deconbenchmark_models.html")

# Get all supported methods
all_methods <- getSupportedMethods()
cat(paste("Found", length(all_methods), "supported methods\n"))

# Initialize data frame for results
model_info <- data.frame(
  Model = character(),
  RequiredInputs = character(),
  ModelType = character(),
  ExpectedOutputs = character(),
  stringsAsFactors = FALSE
)

# Process each method
for (method in all_methods) {
  cat(paste("Processing method:", method, "\n"))
  
  # Get required inputs with error handling
  tryCatch({
    inputs <- getMethodsInputs(method, containerEngine = "singularity")
    required_inputs <- inputs[[method]]
    inputs_str <- paste(required_inputs, collapse = ", ")
    
    # Classify model type with ordered factor for better sorting
    reference_inputs <- c("singleCellExpr", "signature", "cellTypeExpr", "cellTypeExpression")
    marker_inputs <- c("markers")
    basic_inputs <- c("bulk", "nCellTypes")
    
    model_type <- "Unknown"
    if (any(reference_inputs %in% required_inputs)) {
      model_type <- "1-Reference-based"
    } else if (any(marker_inputs %in% required_inputs)) {
      model_type <- "2-Semi-reference-free"
    } else if (all(required_inputs %in% basic_inputs)) {
      model_type <- "3-Reference-free"
    }
    
    # Determine expected outputs
    expected_outputs <- "Cell type proportions (P)"
    if (model_type != "Unknown") {
      expected_outputs <- paste(expected_outputs, ", Signature matrix (S)")
    }
    
    # Add to data frame
    model_info <- rbind(model_info, data.frame(
      Model = method,
      RequiredInputs = inputs_str,
      ModelType = model_type,
      ExpectedOutputs = expected_outputs,
      stringsAsFactors = FALSE
    ))
  }, error = function(e) {
    cat(paste("Error processing method", method, ":", e$message, "\n"))
    model_info <- rbind(model_info, data.frame(
      Model = method,
      RequiredInputs = "Error retrieving inputs",
      ModelType = "4-Unknown",
      ExpectedOutputs = "Unknown",
      stringsAsFactors = FALSE
    ))
  })
}

# Sort by model type (using the prefix numbers for order)
model_info <- model_info[order(model_info$ModelType, model_info$Model), ]

# Clean up display type names by removing the sorting prefix
model_info$DisplayType <- gsub("^[0-9]-", "", model_info$ModelType)

# Count by type for summary
ref_based_count <- sum(grepl("Reference-based", model_info$DisplayType))
ref_free_count <- sum(grepl("Reference-free", model_info$DisplayType))
semi_ref_free_count <- sum(grepl("Semi-reference-free", model_info$DisplayType))

# Generate HTML file
html_content <- paste0('
<!DOCTYPE html>
<html>
<head>
  <title>DeconBenchmark Models Overview</title>
  <style>
    body { font-family: Arial, sans-serif; margin: 20px; }
    table { border-collapse: collapse; width: 100%; margin-top: 20px; }
    th, td { border: 1px solid #ddd; padding: 8px; }
    th { background-color: #f2f2f2; font-weight: bold; position: sticky; top: 0; }
    tr:nth-child(even) { background-color: #f9f9f9; }
    .reference-based { background-color: #e6f7ff; }
    .reference-free { background-color: #f0fff0; }
    .semi-reference-free { background-color: #fff0f0; }
    .unknown { background-color: #f0f0f0; }
    .summary { margin-bottom: 20px; }
    .model-group { background-color: #333; color: white; font-weight: bold; }
  </style>
</head>
<body>
  <h1>DeconBenchmark Models Overview</h1>
  <p>Generated on: ', format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '</p>
  
  <div class="summary">
    <h2>Summary</h2>
    <p>Total models: ', nrow(model_info), '</p>
    <p>Reference-based: ', ref_based_count, 
    ' | Reference-free: ', ref_free_count, 
    ' | Semi-reference-free: ', semi_ref_free_count, '</p>
  </div>
  
  <table>
    <tr>
      <th>Model</th>
      <th>Required Inputs</th>
      <th>Model Type</th>
      <th>Expected Outputs</th>
    </tr>
')

# Track current group to add headers
current_group <- ""

# Add table rows with group headers
for (i in 1:nrow(model_info)) {
  if (model_info$DisplayType[i] != current_group) {
    # Add group header row
    current_group <- model_info$DisplayType[i]
    html_content <- paste0(html_content, '
    <tr class="model-group">
      <td colspan="4">', current_group, ' Models</td>
    </tr>
    ')
  }
  
  # Add model row with appropriate class
  row_class <- tolower(gsub("-", "-", model_info$DisplayType[i]))
  row_class <- gsub(" ", "-", row_class)
  
  html_content <- paste0(html_content, '
    <tr class="', row_class, '">
      <td>', model_info$Model[i], '</td>
      <td>', model_info$RequiredInputs[i], '</td>
      <td>', model_info$DisplayType[i], '</td>
      <td>', model_info$ExpectedOutputs[i], '</td>
    </tr>
  ')
}

# Close HTML
html_content <- paste0(html_content, '
  </table>
  
  <div style="margin-top: 20px;">
    <h3>Legend</h3>
    <ul>
      <li><span style="background-color: #e6f7ff; padding: 2px 5px;">Reference-based</span>: Requires reference data (single-cell expression, signatures, etc.)</li>
      <li><span style="background-color: #f0fff0; padding: 2px 5px;">Reference-free</span>: Only requires bulk data and number of cell types</li>
      <li><span style="background-color: #fff0f0; padding: 2px 5px;">Semi-reference-free</span>: Requires only marker genes</li>
    </ul>
  </div>
</body>
</html>
')

# Write to file
write(html_content, output_file)
cat(paste("HTML report generated:", output_file, "\n"))