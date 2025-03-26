#!/bin/sh
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH --mem 48G
#SBATCH --time 2:00:00
#SBATCH --output=/scratch/lorthiois/logs/%A.o
#SBATCH --error=/scratch/lorthiois/logs/%A.e
#SBATCH --job-name=singleCellExpr_gen

set -e
set -x

# Default parameters (can be overridden by command line arguments)
RLIBRARY="${1:-/work/gr-fe/R_4.3.1}"
SCRIPT="/work/gr-fe/lorthiois/DeconBenchmark/scripts/data/singleCellExpr_generation.R"
INPUT_DATA="${2:-/work/gr-fe/lorthiois/DeconBenchmark/data/pbmc_reference.rds}"
OUTPUT_BASE_DIR="${3:-/work/gr-fe/lorthiois/DeconBenchmark/generated_data}"
# Add the missing mapping file parameter with a default value
MAPPING_FILE="${4:-/work/gr-fe/lorthiois/DeconBenchmark/data/mart_export.txt}"
PREFIX="PBMC"

# Create prefix-specific subdirectory
OUTPUT_DIR="$OUTPUT_BASE_DIR/$PREFIX"
mkdir -p $OUTPUT_DIR

module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/
module load r

start=`date +%s`
echo "START AT $(date)"

Rscript ${SCRIPT} ${RLIBRARY} ${INPUT_DATA} ${OUTPUT_DIR} ${MAPPING_FILE} ${PREFIX}

end=`date +%s`
runtime=$((end-start))
echo "Runtime: $runtime seconds"
echo "Data saved to: $OUTPUT_DIR"