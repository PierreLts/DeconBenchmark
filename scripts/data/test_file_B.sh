#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=0:10:00
#SBATCH --output=/scratch/lorthiois/logs/%A.o
#SBATCH --error=/scratch/lorthiois/logs/%A.e
#SBATCH --job-name=test_read_B

# Exit on error
set -e

# Load R module
module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/
module load r

# Print versions and paths for troubleshooting
echo "=== Environment Information ==="
echo "R version: $(R --version | head -1)"
echo "R library path: $(R -e '.libPaths()' | grep -v '>' | head -1)"
echo "Working directory: $(pwd)"
echo 

# Create temporary R script
TEMP_SCRIPT=$(mktemp)
cat > "$TEMP_SCRIPT" << 'EOF'
# Print system info
cat("R version:", R.version.string, "\n")
cat("Working directory:", getwd(), "\n\n")

# File to test
file_path <- "/work/gr-fe/lorthiois/DeconBenchmark/generated_data/TB1/TB1_singleCellExpr_B.rda"

# Check if file exists
cat("Checking file:", file_path, "\n")
if (!file.exists(file_path)) {
  cat("ERROR: File does not exist!\n")
  quit(status = 1)
}

# Get file info
file_info <- file.info(file_path)
cat("File size:", file_info$size, "bytes\n")
cat("File permissions:", file.access(file_path, 4), "(0=readable)\n")
cat("Last modified:", file_info$mtime, "\n\n")

# Try to read the file header
cat("Testing file header...\n")
tryCatch({
  con <- file(file_path, "rb")
  header <- readBin(con, "raw", 100)
  close(con)
  cat("Successfully read file header\n\n")
}, error = function(e) {
  cat("ERROR reading file header:", e$message, "\n\n")
})

# Try to load the full file
cat("Attempting to load the full file...\n")
tryCatch({
  load_result <- load(file_path)
  cat("File loaded successfully\n")
  cat("Objects in file:", paste(load_result, collapse=", "), "\n")
  
  if (exists("singleCellExpr")) {
    cat("singleCellExpr object found\n")
    cat("Dimensions:", paste(dim(singleCellExpr), collapse=" x "), "\n")
    cat("Class:", class(singleCellExpr), "\n")
    cat("First few row names:", paste(head(rownames(singleCellExpr)), collapse=", "), "\n")
    cat("First few column names:", paste(head(colnames(singleCellExpr)), collapse=", "), "\n")
  } else {
    cat("WARNING: singleCellExpr object not found in file\n")
  }
}, error = function(e) {
  cat("ERROR loading file:", e$message, "\n")
})

# Try to load the A and AB files for comparison
cat("\n=== Comparing with A and AB files ===\n")
for (filter in c("A", "AB")) {
  compare_file <- gsub("_B.rda", paste0("_", filter, ".rda"), file_path)
  cat("Checking", filter, "file:", compare_file, "\n")
  
  if (!file.exists(compare_file)) {
    cat("File does not exist, skipping\n")
    next
  }
  
  file_info <- file.info(compare_file)
  cat("File size:", file_info$size, "bytes\n")
  
  tryCatch({
    load(compare_file)
    if (exists("singleCellExpr")) {
      cat("Successfully loaded singleCellExpr\n")
      cat("Dimensions:", paste(dim(singleCellExpr), collapse=" x "), "\n")
    } else {
      cat("WARNING: singleCellExpr object not found\n")
    }
  }, error = function(e) {
    cat("ERROR loading file:", e$message, "\n")
  })
  cat("\n")
}
EOF

# Run the R script
echo "=== Running R diagnostics ==="
Rscript "$TEMP_SCRIPT"

# Clean up
rm "$TEMP_SCRIPT"

echo "=== Test completed ==="