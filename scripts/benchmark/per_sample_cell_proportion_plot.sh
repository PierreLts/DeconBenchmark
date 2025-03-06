#!/bin/sh
## This script compares predicted proportions with paired sample-specific ground truth
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 32G
#SBATCH --time 1:00:00
#SBATCH --output=/scratch/lorthiois/logs/%A.o
#SBATCH --error=/scratch/lorthiois/logs/%A.e
#SBATCH --job-name=paired_gt_plot

### Useful for debugging
set -e # Exit the script if any statement returns a non-true return value
set -x # Print each line of code being computed

RLIBRARY="/work/gr-fe/R_4.3.1" # R library path
SCRIPT=/work/gr-fe/lorthiois/DeconBenchmark/scripts/benchmark/per_sample_cell_proportion_plot.R
DATASET_PREFIX="${1:-TB}"
METHOD="${2:-CIBERSORT}"
OUTPUT_BASE_DIR="/work/gr-fe/lorthiois/DeconBenchmark/benchmark_results"
INCLUDE_OVERALL_GT="${3:-TRUE}"

module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/
module load r

start=`date +%s`
echo "START AT $(date)"

# Create output directory
OUTPUT_DIR="${OUTPUT_BASE_DIR}/${DATASET_PREFIX}/${METHOD}"
mkdir -p ${OUTPUT_DIR}

# Run the R script (ORDER IS IMPORTANT)
Rscript ${SCRIPT} ${RLIBRARY} ${DATASET_PREFIX} ${METHOD} ${OUTPUT_BASE_DIR} ${INCLUDE_OVERALL_GT}

# Print end date and echo total runtime
end=`date +%s`
runtime=$((end-start))
echo "Runtime: $runtime seconds"
echo "Benchmarking results saved to: ${OUTPUT_DIR}"