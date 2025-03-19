#!/bin/sh
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 16G
#SBATCH --time 0:30:00
#SBATCH --output=/scratch/lorthiois/logs/%A.o
#SBATCH --error=/scratch/lorthiois/logs/%A.e
#SBATCH --job-name=BE_GT_sample_gen
### WARNING: THE FILE IS INCORRECT BUT OUTPUTS A FORMMAT_VALID FILE FOR BENCHMARKING
# Enhanced error handling
set -e  # Exit immediately if a command exits with a non-zero status
set -x  # Print commands and their arguments as they are executed

# Default parameters (can be overridden by command line arguments)
RLIBRARY="${1:-/work/gr-fe/R_4.3.1}"
OUTPUT_BASE_DIR="${2:-/work/gr-fe/lorthiois/DeconBenchmark/generated_data}"
PREFIX="${3:-BE}"  # Default prefix is BE, but can be overridden

# Create prefix-specific subdirectory
OUTPUT_DIR="$OUTPUT_BASE_DIR/$PREFIX"
echo "Output directory: $OUTPUT_DIR"

# Check if directory exists, create if it doesn't
if [ ! -d "$OUTPUT_DIR" ]; then
    echo "Creating output directory: $OUTPUT_DIR"
    mkdir -p "$OUTPUT_DIR"
    if [ $? -ne 0 ]; then
        echo "ERROR: Failed to create output directory"
        exit 1
    fi
else
    echo "Output directory already exists"
    # Test write permission
    touch "$OUTPUT_DIR/test_write_permission" 2>/dev/null
    if [ $? -ne 0 ]; then
        echo "ERROR: No write permission in output directory"
        exit 1
    else
        echo "Write permission confirmed"
        rm "$OUTPUT_DIR/test_write_permission"
    fi
fi

# Path to R script
SCRIPT="/work/gr-fe/lorthiois/DeconBenchmark/scripts/data/BE_GT_per_sample_generation.R"

# Verify script exists
if [ ! -f "$SCRIPT" ]; then
    echo "ERROR: R script not found at: $SCRIPT"
    exit 1
else
    echo "R script found at: $SCRIPT"
fi

# Load R module
echo "Loading R module..."
module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/
module load r

# Test R installation
Rscript --version
if [ $? -ne 0 ]; then
    echo "ERROR: R is not available or not working properly"
    exit 1
fi

start=`date +%s`
echo "START AT $(date)"
echo "Using prefix: $PREFIX"

# Run R script
echo "Running R script: $SCRIPT"
Rscript ${SCRIPT} ${RLIBRARY} ${OUTPUT_DIR} ${PREFIX}
SCRIPT_EXIT_CODE=$?

if [ $SCRIPT_EXIT_CODE -ne 0 ]; then
    echo "ERROR: R script execution failed with exit code $SCRIPT_EXIT_CODE"
    exit $SCRIPT_EXIT_CODE
fi

# Verify output files were created
RDA_FILE="$OUTPUT_DIR/${PREFIX}_GT_proportions_per_sample.rda"
CSV_FILE="$OUTPUT_DIR/${PREFIX}_GT_proportions_per_sample.csv"

if [ ! -f "$RDA_FILE" ]; then
    echo "WARNING: Expected RDA file not found: $RDA_FILE"
else
    echo "RDA file created: $RDA_FILE ($(stat -c%s "$RDA_FILE") bytes)"
fi

if [ ! -f "$CSV_FILE" ]; then
    echo "WARNING: Expected CSV file not found: $CSV_FILE"
else
    echo "CSV file created: $CSV_FILE ($(stat -c%s "$CSV_FILE") bytes)"
fi

end=`date +%s`
runtime=$((end-start))
echo "Runtime: $runtime seconds"
echo "Per-sample ground truth data processing completed for prefix: $PREFIX"