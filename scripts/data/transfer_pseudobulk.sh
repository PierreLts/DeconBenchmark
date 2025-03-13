#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=0:30:00
#SBATCH --job-name=pseudobulk_transfer

set -e  # Exit on error
set -x  # Print commands as they're executed

# Default parameters (can be overridden by command line arguments)
RLIBRARY="${1:-/work/gr-fe/R_4.3.1}"
SCRIPT="/work/gr-fe/lorthiois/DeconBenchmark/scripts/data/transfer_pseudobulk.R"
SOURCE_FILE="${2:-/work/gr-fe/lorthiois/DeconBenchmark/data/pseudobulk_counts_120k.csv}"
MAPPING_FILE="${3:-/work/gr-fe/lorthiois/DeconBenchmark/data/mart_export.txt}"
OUTPUT_BASE_DIR="${4:-/work/gr-fe/lorthiois/DeconBenchmark/generated_data}"
PREFIX="${5:-TB1}"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_BASE_DIR/$PREFIX"

# Load R module
module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/
module load r

start=`date +%s`
echo "START AT $(date)"
echo "Processing pseudobulk data for $PREFIX"

# Run the R script (removed the sample filter parameter)
Rscript ${SCRIPT} ${RLIBRARY} ${SOURCE_FILE} ${MAPPING_FILE} ${OUTPUT_BASE_DIR} ${PREFIX}

end=`date +%s`
runtime=$((end-start))
echo "Runtime: $runtime seconds"
echo "Pseudobulk data processed and saved to: $OUTPUT_BASE_DIR/$PREFIX"