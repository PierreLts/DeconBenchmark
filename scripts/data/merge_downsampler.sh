#!/bin/sh
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH --mem 64G
#SBATCH --time 2:00:00
#SBATCH --output=/scratch/lorthiois/logs/%A.o
#SBATCH --error=/scratch/lorthiois/logs/%A.e
#SBATCH --job-name=merge_downsample

set -e  # Exit on error
set -x  # Print commands as they're executed

# ===== PARAMETERS TO ADJUST =====
OUTPUT_PREFIX="TB"
CELLS_PER_TYPE=25 # 4 batches --> multiply by 4 nb of cells/celltype
FILTER_TYPE="AB"  # Default filter type (AB, A, or B)
BATCH1="TB1"
BATCH2="TB2"
BATCH3="TB3"
BATCH4="TB4"
# ================================

# Default parameters
RLIBRARY="${1:-/work/gr-fe/R_4.3.1}"
SCRIPT_DIR="/work/gr-fe/lorthiois/DeconBenchmark/scripts/data"
INPUT_BASE_DIR="${2:-/work/gr-fe/lorthiois/DeconBenchmark/generated_data}"
OUTPUT_BASE_DIR="${3:-/work/gr-fe/lorthiois/DeconBenchmark/generated_data}"

# Override parameters if provided as command line arguments
if [ -n "$4" ]; then
    OUTPUT_PREFIX="$4"
fi
if [ -n "$5" ]; then
    CELLS_PER_TYPE="$5"
fi
if [ -n "$6" ]; then
    FILTER_TYPE="$6"
fi
if [ -n "$7" ]; then
    BATCH1="$7"
fi
if [ -n "$8" ]; then
    BATCH2="$8"
fi
if [ -n "$9" ]; then
    BATCH3="$9"
fi
if [ -n "${10}" ]; then
    BATCH4="${10}"
fi

# Validate filter type
if [ "$FILTER_TYPE" != "A" ] && [ "$FILTER_TYPE" != "B" ] && [ "$FILTER_TYPE" != "AB" ]; then
    echo "ERROR: Filter type must be 'A', 'B', or 'AB'. Got '$FILTER_TYPE'"
    exit 1
fi

# Calculate number of batches and total cells
NUM_BATCHES=0
BATCHES=""

if [ -n "$BATCH1" ]; then
    NUM_BATCHES=$((NUM_BATCHES + 1))
    BATCHES="$BATCHES $BATCH1"
fi
if [ -n "$BATCH2" ]; then
    NUM_BATCHES=$((NUM_BATCHES + 1))
    BATCHES="$BATCHES $BATCH2"
fi
if [ -n "$BATCH3" ]; then
    NUM_BATCHES=$((NUM_BATCHES + 1))
    BATCHES="$BATCHES $BATCH3"
fi
if [ -n "$BATCH4" ]; then
    NUM_BATCHES=$((NUM_BATCHES + 1))
    BATCHES="$BATCHES $BATCH4"
fi

TOTAL_CELLS=$((NUM_BATCHES * CELLS_PER_TYPE))

# Set output info
OUTPUT_PREFIX_WITH_SUFFIX="${OUTPUT_PREFIX}_D${TOTAL_CELLS}"
OUTPUT_DIR="${OUTPUT_BASE_DIR}/${OUTPUT_PREFIX_WITH_SUFFIX}"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Load R module
module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/
module load r

echo "=== Merging and Downsampling Datasets ==="
echo "Output prefix: ${OUTPUT_PREFIX} -> ${OUTPUT_PREFIX_WITH_SUFFIX}"
echo "Cells per type per batch: ${CELLS_PER_TYPE}"
echo "Total cells per type: ${TOTAL_CELLS}"
echo "Number of batches: ${NUM_BATCHES}"
echo "Batches: ${BATCHES}"
echo "Filter type: ${FILTER_TYPE}"

# Run the R script
start=`date +%s`
echo "START AT $(date)"

R_SCRIPT="${SCRIPT_DIR}/merge_downsampler.R"

# Check if the script exists
if [ ! -f "$R_SCRIPT" ]; then
    echo "ERROR: R script not found at $R_SCRIPT"
    exit 1
fi

Rscript ${R_SCRIPT} ${RLIBRARY} ${OUTPUT_PREFIX} ${CELLS_PER_TYPE} ${FILTER_TYPE} ${BATCH1} ${BATCH2} ${BATCH3} ${BATCH4}
RSCRIPT_STATUS=$?

if [ $RSCRIPT_STATUS -ne 0 ]; then
    echo "ERROR: R script failed with status $RSCRIPT_STATUS"
    exit $RSCRIPT_STATUS
fi

end=`date +%s`
runtime=$((end-start))
echo "Runtime: $runtime seconds"
echo "Merged and downsampled dataset created: ${OUTPUT_PREFIX_WITH_SUFFIX}"
echo "Files saved to: ${OUTPUT_DIR}"