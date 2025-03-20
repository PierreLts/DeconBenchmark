#!/bin/sh
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 32G
#SBATCH --time 2:00:00
#SBATCH --output=/scratch/lorthiois/logs/%A.o
#SBATCH --error=/scratch/lorthiois/logs/%A.e
#SBATCH --job-name=downsampler

set -e  # Exit on error
set -x  # Print commands as they're executed

# ===== PARAMETERS TO ADJUST =====
PREFIX="TB1"               # Input dataset prefix
CELLS_PER_TYPE=100           # Number of cells to keep per cell type
# ================================

# Default parameters
RLIBRARY="${1:-/work/gr-fe/R_4.3.1}"
SCRIPT_DIR="/work/gr-fe/lorthiois/DeconBenchmark/scripts/data"
INPUT_BASE_DIR="${2:-/work/gr-fe/lorthiois/DeconBenchmark/generated_data}"
OUTPUT_BASE_DIR="${3:-/work/gr-fe/lorthiois/DeconBenchmark/generated_data}"

# Override PREFIX if provided as 4th argument
if [ -n "$4" ]; then
    PREFIX="$4"
fi
# Override CELLS_PER_TYPE if provided as 5th argument
if [ -n "$5" ]; then
    CELLS_PER_TYPE="$5"
fi

# Set paths
INPUT_DIR="${INPUT_BASE_DIR}/${PREFIX}"
OUTPUT_PREFIX="${PREFIX}_D${CELLS_PER_TYPE}"
OUTPUT_DIR="${OUTPUT_BASE_DIR}/${OUTPUT_PREFIX}"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Load R module
module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/
module load r

echo "=== Downsampling Dataset ==="
echo "Source: ${INPUT_DIR}"
echo "Destination: ${OUTPUT_DIR}" 
echo "Prefix: ${PREFIX} -> ${OUTPUT_PREFIX}"
echo "Cells per type: ${CELLS_PER_TYPE}"

# First copy all files from source to destination (except the ones we'll downsample)
echo "Copying all files..."
for file in "${INPUT_DIR}"/*; do
    if [ -f "$file" ]; then
        filename=$(basename "$file")
        
        # Skip the files we'll downsample - we'll handle these separately
        if [[ "$filename" == *"_singleCellExpr"* ]] || [[ "$filename" == *"_singleCellLabels"* ]] || [[ "$filename" == *"_singleCellSubjects"* ]]; then
            echo "Skipping for downsampling: $filename"
            continue
        fi
        
        # Rename file with new prefix and copy
        new_filename="${filename/$PREFIX/$OUTPUT_PREFIX}"
        cp "$file" "${OUTPUT_DIR}/${new_filename}"
        echo "Copied: $filename -> ${new_filename}"
    fi
done

# Run the R script to downsample the single cell data
echo "Running downsampling..."
R_SCRIPT="${SCRIPT_DIR}/downsampler.R"

Rscript ${R_SCRIPT} ${RLIBRARY} ${INPUT_DIR} ${OUTPUT_DIR} ${PREFIX} ${CELLS_PER_TYPE}
RSCRIPT_STATUS=$?

if [ $RSCRIPT_STATUS -ne 0 ]; then
    echo "ERROR: R script failed with status $RSCRIPT_STATUS"
    exit $RSCRIPT_STATUS
fi

echo "Renaming downsampled files..."
# Rename all files created by the R script from PREFIX to OUTPUT_PREFIX
for ext in "rda" "csv"; do
    for file_pattern in "_singleCellExpr" "_singleCellLabels" "_singleCellSubjects"; do
        # Find files matching the pattern with both extensions
        for file in $(find "${OUTPUT_DIR}" -name "${PREFIX}${file_pattern}*${ext}" -type f); do
            if [ -f "$file" ]; then
                filename=$(basename "$file")
                
                # Rename file from PREFIX to OUTPUT_PREFIX
                new_filename="${filename/$PREFIX/$OUTPUT_PREFIX}"
                mv "$file" "${OUTPUT_DIR}/${new_filename}"
                echo "Renamed: $filename -> ${new_filename}"
            fi
        done
    done
done

echo "Downsampling completed!"
echo "New dataset created: ${OUTPUT_PREFIX}"