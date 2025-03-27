#!/bin/sh
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 32G
#SBATCH --time 1:00:00
#SBATCH --output=/scratch/lorthiois/logs/%A.o
#SBATCH --error=/scratch/lorthiois/logs/%A.e
#SBATCH --job-name=common_genes

set -e  # Exit on error
set -x  # Print commands as they're executed

# Default parameters
RLIBRARY="${1:-/work/gr-fe/R_4.3.1}"
INPUT_FILE1="${2:-/work/gr-fe/lorthiois/DeconBenchmark/generated_data/TB_D100/TB_D100_singleCellExpr_AB.rda}"
INPUT_FILE2="${3:-/work/gr-fe/lorthiois/DeconBenchmark/generated_data/PBMC_D100/PBMC_D100_singleCellExpr_AB.rda}"
BASE_DIR="${4:-/work/gr-fe/lorthiois/DeconBenchmark/generated_data}"

# Handle relative or partial paths
if [[ "$INPUT_FILE1" != /* ]]; then
  # If it's not an absolute path, prepend the base directory
  INPUT_FILE1="${BASE_DIR}/${INPUT_FILE1}"
fi
if [[ "$INPUT_FILE2" != /* ]]; then
  # If it's not an absolute path, prepend the base directory
  INPUT_FILE2="${BASE_DIR}/${INPUT_FILE2}"
fi

# Check if input files exist
if [ ! -f "$INPUT_FILE1" ]; then
  echo "ERROR: Input file 1 not found: $INPUT_FILE1"
  exit 1
fi
if [ ! -f "$INPUT_FILE2" ]; then
  echo "ERROR: Input file 2 not found: $INPUT_FILE2"
  exit 1
fi

# Extract file information
INPUT1_DIR=$(dirname "$INPUT_FILE1")
INPUT2_DIR=$(dirname "$INPUT_FILE2")
INPUT1_PREFIX=$(basename "$INPUT1_DIR")
INPUT2_PREFIX=$(basename "$INPUT2_DIR")
INPUT1_BASENAME=$(basename "$INPUT_FILE1")
INPUT2_BASENAME=$(basename "$INPUT_FILE2")

# Create output paths
OUTPUT1_PREFIX="${INPUT1_PREFIX}-common"
OUTPUT2_PREFIX="${INPUT2_PREFIX}-common"
# Make sure we only replace the prefix at the beginning of the filename
OUTPUT1_BASENAME=$(echo "$INPUT1_BASENAME" | sed "s/^${INPUT1_PREFIX}/${OUTPUT1_PREFIX}/")
OUTPUT2_BASENAME=$(echo "$INPUT2_BASENAME" | sed "s/^${INPUT2_PREFIX}/${OUTPUT2_PREFIX}/")
OUTPUT1_DIR="${BASE_DIR}/${OUTPUT1_PREFIX}"
OUTPUT2_DIR="${BASE_DIR}/${OUTPUT2_PREFIX}"
OUTPUT_FILE1="${OUTPUT1_DIR}/${OUTPUT1_BASENAME}"
OUTPUT_FILE2="${OUTPUT2_DIR}/${OUTPUT2_BASENAME}"

# Create output directories
mkdir -p "$OUTPUT1_DIR"
mkdir -p "$OUTPUT2_DIR"

# Print summary information
echo "=== Summary ==="
echo "Input files processed:"
echo "  1. $INPUT_FILE1"
echo "  2. $INPUT_FILE2"
echo ""
echo "Output files created:"
echo "  1. $OUTPUT_FILE1 (RDA)"
echo "  2. ${OUTPUT_FILE1%.*}_50x50.csv (CSV preview)"
echo "  3. $OUTPUT_FILE2 (RDA)"
echo "  4. ${OUTPUT_FILE2%.*}_50x50.csv (CSV preview)"

# Path to the R script - it should be in the same directory as this script
# Specify the absolute path to the R script
R_SCRIPT="/work/gr-fe/lorthiois/DeconBenchmark/scripts/data/common_genes.R"

# Print R script path for debugging
echo "R script path: $R_SCRIPT"

# Verify that the R script exists
if [ ! -f "$R_SCRIPT" ]; then
  echo "ERROR: R script not found at: $R_SCRIPT"
  echo "Please ensure the R script is placed at the correct location."
  exit 1
fi

# Check if R script exists
if [ ! -f "$R_SCRIPT" ]; then
  echo "ERROR: R script not found: $R_SCRIPT"
  exit 1
fi

# Load R module
module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/
module load r

start=`date +%s`
echo "START AT $(date)"

# Run the R script
Rscript "${R_SCRIPT}" "${RLIBRARY}" "${INPUT_FILE1}" "${INPUT_FILE2}" "${OUTPUT_FILE1}" "${OUTPUT_FILE2}"

# Check exit status
if [ $? -ne 0 ]; then
  echo "ERROR: R script execution failed"
  exit 1
fi

end=`date +%s`
runtime=$((end-start))
echo "Runtime: $runtime seconds"
echo "Processing complete!"