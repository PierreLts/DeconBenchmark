#!/bin/sh
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 16G
#SBATCH --time 0:30:00
#SBATCH --output=/scratch/lorthiois/logs/%A.o
#SBATCH --error=/scratch/lorthiois/logs/%A.e
#SBATCH --job-name=ground_truth

set -e
set -x

# Default parameters (can be overridden by command line arguments)
RLIBRARY="${1:-/work/gr-fe/R_4.3.1}"
SCRIPT="/work/gr-fe/lorthiois/DeconBenchmark/scripts/data/ground_truth.R"
INPUT_DIR="${2:-/work/gr-fe/lorthiois/DeconBenchmark/generated_data}"
OUTPUT_BASE_DIR="${3:-/work/gr-fe/lorthiois/DeconBenchmark/generated_data}"
PREFIX="${4:-TB}"

# Create prefix-specific subdirectory
OUTPUT_DIR="$OUTPUT_BASE_DIR/$PREFIX"
mkdir -p $OUTPUT_DIR

module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/
module load r

start=`date +%s`
echo "START AT $(date)"

Rscript ${SCRIPT} ${RLIBRARY} ${INPUT_DIR} ${OUTPUT_DIR} ${PREFIX}

end=`date +%s`
runtime=$((end-start))
echo "Runtime: $runtime seconds"
echo "Data saved to: $OUTPUT_DIR"