#!/bin/sh
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH --mem 48G
#SBATCH --time 3:00:00
#SBATCH --output=/scratch/lorthiois/logs/%A.o
#SBATCH --error=/scratch/lorthiois/logs/%A.e
#SBATCH --job-name=singleCellExpr_AB_gen

set -e
set -x

# ===== EXPLICIT SETTINGS - MODIFY THESE DIRECTLY =====
# Output prefix for generated files
PREFIX="PBMC"
# Sample filter: "A" for only A samples, "B" for only B samples, "AB" for both
SAMPLE_FILTER="AB"
# =====

# Default parameters (can be overridden by command line arguments)
RLIBRARY="${1:-/work/gr-fe/R_4.3.1}"
SCRIPT="/work/gr-fe/lorthiois/DeconBenchmark/scripts/data/singleCellExpr_AB_generation.R"
INPUT_DATA="${2:-/work/gr-fe/lorthiois/DeconBenchmark/data/pbmc_reference.rds}"
OUTPUT_BASE_DIR="${3:-/work/gr-fe/lorthiois/DeconBenchmark/generated_data}"
MAPPING_FILE="${4:-/work/gr-fe/lorthiois/DeconBenchmark/data/mart_export.txt}"

# Override defaults with command line arguments if provided
if [ "$5" != "" ]; then
    PREFIX="$5"
fi
if [ "$6" != "" ]; then
    SAMPLE_FILTER="$6"
fi

# Create prefix-specific subdirectory
OUTPUT_DIR="$OUTPUT_BASE_DIR/$PREFIX"
mkdir -p $OUTPUT_DIR

# Log settings
echo "===== SETTINGS ====="
echo "PREFIX: $PREFIX"
echo "SAMPLE_FILTER: $SAMPLE_FILTER"
echo "OUTPUT_DIR: $OUTPUT_DIR"
echo "INPUT_DATA: $INPUT_DATA"
echo "===================="

module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/
module load r

start=`date +%s`
echo "START AT $(date)"
echo "Filtering samples by: $SAMPLE_FILTER"

Rscript ${SCRIPT} ${RLIBRARY} ${INPUT_DATA} ${OUTPUT_DIR} ${MAPPING_FILE} ${PREFIX} ${SAMPLE_FILTER}

end=`date +%s`
runtime=$((end-start))
echo "Runtime: $runtime seconds"
echo "Data saved to: $OUTPUT_DIR with prefix $PREFIX and filter $SAMPLE_FILTER"