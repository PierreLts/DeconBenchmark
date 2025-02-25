#!/bin/sh
## This script compares performance metrics across all deconvolution models
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 2
#SBATCH --mem 8G
#SBATCH --time 1:00:00
#SBATCH --output=/scratch/lorthiois/logs/%A.o
#SBATCH --error=/scratch/lorthiois/logs/%A.e
#SBATCH --job-name=compare_models

### Useful for debugging
set -e # Exit the script if any statement returns a non-true return value
set -x # Print each line of code being computed

RLIBRARY="/work/gr-fe/R_4.3.1" # R library path
SCRIPT=/work/gr-fe/lorthiois/DeconBenchmark/scripts/compare_models.R
BENCHMARK_DIR="/work/gr-fe/lorthiois/DeconBenchmark/benchmark_results"
OUTPUT_DIR="/work/gr-fe/lorthiois/DeconBenchmark/comparison_results"

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/
module load r

start=`date +%s`
echo "START AT $(date)"

# Run the R script (ORDER IS IMPORTANT)
Rscript ${SCRIPT} ${RLIBRARY} ${BENCHMARK_DIR} ${OUTPUT_DIR}

# Print end date and echo total runtime
end=`date +%s`
runtime=$((end-start))
echo "Runtime: $runtime seconds"
echo "Comparison results saved to: ${OUTPUT_DIR}"