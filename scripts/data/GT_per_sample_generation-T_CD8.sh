#!/bin/sh
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 32G
#SBATCH --time 1:00:00
#SBATCH --output=/scratch/lorthiois/logs/%A.o
#SBATCH --error=/scratch/lorthiois/logs/%A.e
#SBATCH --job-name=GT_per_sample_T_CD8_gen

set -e  # Exit on error
set -x  # Print commands as they're executed

# Default parameters (can be overridden by command line arguments)
RLIBRARY="${1:-/work/gr-fe/R_4.3.1}"
SCRIPT="/work/gr-fe/lorthiois/DeconBenchmark/scripts/data/GT_per_sample_generation-T_CD8.R"
INPUT_DATA="${2:-/work/gr-fe/lorthiois/DeconBenchmark/data/merged_batches.rds}"
OUTPUT_BASE_DIR="${3:-/work/gr-fe/lorthiois/DeconBenchmark/generated_data}"
PREFIX="PBMC"

# Create prefix-specific subdirectory
OUTPUT_DIR="$OUTPUT_BASE_DIR/$PREFIX"
mkdir -p $OUTPUT_DIR

module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/
module load r

start=`date +%s`
echo "START AT $(date)"
echo "Generating per-sample ground truth with combined T_CD8 cells"

Rscript ${SCRIPT} ${RLIBRARY} ${INPUT_DATA} ${OUTPUT_DIR} ${PREFIX}

end=`date +%s`
runtime=$((end-start))
echo "Runtime: $runtime seconds"
echo "Data saved to: $OUTPUT_DIR"