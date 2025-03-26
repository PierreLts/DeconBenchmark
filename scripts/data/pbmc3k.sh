#!/bin/sh
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 32G
#SBATCH --time 1:00:00
#SBATCH --output=/scratch/lorthiois/logs/%A.o
#SBATCH --error=/scratch/lorthiois/logs/%A.e
#SBATCH --job-name=pbmc3k

set -e  # Exit on error
set -x  # Print commands as they're executed

# Default parameters
RLIBRARY="${1:-/work/gr-fe/R_4.3.1}"
SCRIPT_PATH="${2:-/work/gr-fe/lorthiois/DeconBenchmark/scripts/data/pbmc3k.R}"
INPUT_DATA="${3:-/work/gr-fe/lorthiois/DeconBenchmark/data/hg19/}"
OUTPUT_FILE="${4:-/work/gr-fe/lorthiois/DeconBenchmark/data/pbmc_reference.rds}"

# Print usage if help flag is provided
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    echo "Usage: $0 [R_LIBRARY_PATH] [SCRIPT_PATH] [INPUT_DATA_PATH] [OUTPUT_FILE_PATH]"
    echo ""
    echo "Parameters:"
    echo "  R_LIBRARY_PATH    - Path to R library (default: /work/gr-fe/R_4.3.1)"
    echo "  SCRIPT_PATH       - Path to pbmc3k.R script (default: /work/gr-fe/lorthiois/DeconBenchmark/scripts/data/pbmc3k.R)"
    echo "  INPUT_DATA_PATH   - Path to input 10X data directory (default: /work/gr-fe/lorthiois/DeconBenchmark/data/hg19/)"
    echo "  OUTPUT_FILE_PATH  - Path to save the output RDS file (default: /work/gr-fe/lorthiois/DeconBenchmark/data/pbmc_reference.rds)"
    exit 0
fi

# Ensure input directory exists
if [ ! -d "$INPUT_DATA" ]; then
    echo "ERROR: Input directory does not exist: $INPUT_DATA"
    exit 1
fi

# Ensure output directory exists
OUTPUT_DIR=$(dirname "$OUTPUT_FILE")
mkdir -p "$OUTPUT_DIR"

# Load R module
module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/
module load r

start=`date +%s`
echo "START AT $(date)"
echo "Processing PBMC dataset from: $INPUT_DATA"
echo "Output will be saved to: $OUTPUT_FILE"

# Run the R script with the provided arguments
Rscript ${SCRIPT_PATH} ${RLIBRARY} ${INPUT_DATA} ${OUTPUT_FILE}

# Check for successful completion
if [ $? -ne 0 ]; then
    echo "ERROR: R script execution failed"
    exit 1
fi

# Print end date and runtime
end=`date +%s`
runtime=$((end-start))
echo "Runtime: $runtime seconds"
echo "PBMC processing completed"
echo "Results saved to: $OUTPUT_FILE"