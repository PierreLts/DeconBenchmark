#!/bin/sh
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 32
#SBATCH --mem 128G
#SBATCH --time 1:00:00
#SBATCH --output=/scratch/lorthiois/logs/%A.o
#SBATCH --error=/scratch/lorthiois/logs/%A.e
#SBATCH --job-name=pseudobulk_gen

set -e  # Exit on error
set -x  # Print commands as they're executed

# ===== EXPLICIT SETTINGS - MODIFY THESE DIRECTLY =====
# Output prefix for generated files
PREFIX="TB1"
# Sample filter: "A" for only A samples, "B" for only B samples, "AB" for both
SAMPLE_FILTER="AB"
# =====

# Default parameters (can be overridden by command line arguments)
RLIBRARY="${1:-/work/gr-fe/R_4.3.1}"
SCRIPT="/work/gr-fe/lorthiois/DeconBenchmark/scripts/data/pseudobulk.R"
SEURAT_FILE="${2:-/work/gr-fe/lorthiois/DeconBenchmark/data/GFB-33245_HFKJMDSXC_2_scRNAseqWTATBseverityrun1_Seurat.rds}"
MAPPING_FILE="${3:-/work/gr-fe/lorthiois/DeconBenchmark/data/mart_export.txt}"
OUTPUT_BASE_DIR="${4:-/work/gr-fe/lorthiois/DeconBenchmark/generated_data}"

# Override defaults with command line arguments if provided
if [ "$5" != "" ]; then
    PREFIX="$5"
fi
if [ "$6" != "" ]; then
    SAMPLE_FILTER="$6"
fi

# Create prefix-specific subdirectory
OUTPUT_DIR="$OUTPUT_BASE_DIR"
mkdir -p $OUTPUT_DIR

# Log settings
echo "===== SETTINGS ====="
echo "PREFIX: $PREFIX"
echo "SAMPLE_FILTER: $SAMPLE_FILTER"
echo "OUTPUT_DIR: $OUTPUT_DIR"
echo "SEURAT_FILE: $SEURAT_FILE"
echo "MAPPING_FILE: $MAPPING_FILE"
echo "===================="

module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/
module load r

start=`date +%s`
echo "START AT $(date)"
echo "Generating pseudobulk data with filter: $SAMPLE_FILTER"

Rscript ${SCRIPT} ${RLIBRARY} ${SEURAT_FILE} ${MAPPING_FILE} ${OUTPUT_DIR} ${PREFIX} ${SAMPLE_FILTER}

end=`date +%s`
runtime=$((end-start))
echo "Runtime: $runtime seconds"
echo "Pseudobulk data saved to: $OUTPUT_DIR/$PREFIX with filter $SAMPLE_FILTER"