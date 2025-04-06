#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=0:30:00
#SBATCH --output=/scratch/lorthiois/logs/%A.o
#SBATCH --error=/scratch/lorthiois/logs/%A.e
#SBATCH --job-name=model_table

set -e  # Exit on error
set -x  # Print commands as they're executed

# Default parameters
RLIBRARY="/work/gr-fe/R_4.3.1"
SCRIPT_DIR="/work/gr-fe/lorthiois/DeconBenchmark/scripts/benchmark"
INPUT_CSV="/work/gr-fe/lorthiois/DeconBenchmark/data/Models.csv"
OUTPUT_FILE="/work/gr-fe/lorthiois/DeconBenchmark/benchmark_results/models_comparison.html"
BENCHMARK_DIR="/work/gr-fe/lorthiois/DeconBenchmark/benchmark_results"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    -i|--input)
      INPUT_CSV="$2"
      shift 2
      ;;
    -o|--output)
      OUTPUT_FILE="$2"
      shift 2
      ;;
    -b|--benchmark-dir)
      BENCHMARK_DIR="$2"
      shift 2
      ;;
    -l|--library)
      RLIBRARY="$2"
      shift 2
      ;;
    -h|--help)
      echo "Usage: $(basename $0) [options]"
      echo "Options:"
      echo "  -i, --input FILE          Input CSV with model requirements (default: $INPUT_CSV)"
      echo "  -o, --output FILE         Output file path (default: $OUTPUT_FILE)"
      echo "  -b, --benchmark-dir DIR   Directory containing benchmark results (default: $BENCHMARK_DIR)"
      echo "  -l, --library PATH        R library path (default: $RLIBRARY)"
      echo "  -h, --help                Show this help message"
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

# Create output directory if needed
OUTPUT_DIR=$(dirname "$OUTPUT_FILE")
mkdir -p "$OUTPUT_DIR"

# Load R module
module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/
module load r

start=`date +%s`
echo "START AT $(date)"
echo "Input CSV: $INPUT_CSV"
echo "Output file: $OUTPUT_FILE"
echo "Benchmark directory: $BENCHMARK_DIR"

# Install required R packages if needed
Rscript - <<EOF
# Check and install required packages
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")
if (!requireNamespace("knitr", quietly = TRUE)) install.packages("knitr")
if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr")
EOF

# Check if input file exists
if [ ! -f "$INPUT_CSV" ]; then
  echo "ERROR: Input CSV file not found: $INPUT_CSV"
  exit 1
fi

# Check if R script exists
if [ ! -f "${SCRIPT_DIR}/table.R" ]; then
  echo "ERROR: R script not found: ${SCRIPT_DIR}/table.R"
  exit 1
fi

# Run the R script
Rscript "${SCRIPT_DIR}/table.R" "$INPUT_CSV" "$OUTPUT_FILE" "$BENCHMARK_DIR"

# Check exit status
if [ $? -ne 0 ]; then
  echo "ERROR: R script execution failed"
  exit 1
fi

# Also generate a simple CSV version for easier data processing
CSV_OUTPUT="${OUTPUT_FILE%.*}.csv"
echo "CSV version available at: $CSV_OUTPUT"

end=`date +%s`
runtime=$((end-start))
echo "Runtime: $runtime seconds"
echo "Table generation completed successfully!"
echo "HTML output: $OUTPUT_FILE"
echo "CSV output: $CSV_OUTPUT"