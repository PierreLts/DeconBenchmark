#!/bin/sh
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH --mem 48G
#SBATCH --time 2:00:00
#SBATCH --output=/scratch/lorthiois/logs/%A.o
#SBATCH --error=/scratch/lorthiois/logs/%A.e
#SBATCH --job-name=pbmc_processor

set -e  # Exit on error
set -x  # Print commands as they're executed

# Default parameters (customize these as needed)
RLIBRARY="/work/gr-fe/R_4.3.1"
INPUT_DATA="/work/gr-fe/lorthiois/DeconBenchmark/data/pbmc_reference.rds"
OUTPUT_DIR="/work/gr-fe/lorthiois/DeconBenchmark/generated_data"
MAPPING_FILE="/work/gr-fe/lorthiois/DeconBenchmark/data/mart_export.txt"
PREFIX="PBMC"

# Override defaults with command line arguments if provided
if [ "$1" != "" ]; then
    RLIBRARY="$1"
fi
if [ "$2" != "" ]; then
    INPUT_DATA="$2"
fi
if [ "$3" != "" ]; then
    OUTPUT_DIR="$3"
fi
if [ "$4" != "" ]; then
    MAPPING_FILE="$4"
fi
if [ "$5" != "" ]; then
    PREFIX="$5"
fi

# Create the destination directory
DEST_DIR="$OUTPUT_DIR/$PREFIX"
mkdir -p "$DEST_DIR"

SCRIPT_PATH="/work/gr-fe/lorthiois/DeconBenchmark/scripts/data/PBMC_singleCellExpr.R"

# Load R module
module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/
module load r

echo "=== Starting PBMC Processing ==="
echo "Input data: $INPUT_DATA"
echo "Output directory: $DEST_DIR"
echo "Prefix: $PREFIX"

start=`date +%s`
echo "START AT $(date)"

# Run the R script
Rscript ${SCRIPT_PATH} ${RLIBRARY} ${INPUT_DATA} ${DEST_DIR} ${MAPPING_FILE} ${PREFIX}

end=`date +%s`
runtime=$((end-start))
echo "Runtime: $runtime seconds"
echo "Data saved to: $DEST_DIR"