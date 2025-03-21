#!/bin/sh
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 32G
#SBATCH --time 00:15:00
#SBATCH --output=/scratch/lorthiois/logs/%A.o
#SBATCH --error=/scratch/lorthiois/logs/%A.e
#SBATCH --job-name=detailed_metrics

RLIBRARY="/work/gr-fe/R_4.3.1" 
SCRIPT=/work/gr-fe/lorthiois/DeconBenchmark/scripts/benchmark/detailed_metrics.R
DATASET_PREFIX="${1:-TB}"
INCLUDE_OVERALL_GT="${2:-TRUE}"
SAMPLE_FILTER="${3:-AB}"  # Parameter for sample filter
GLOBAL_LOG_DIR="${4:-/work/gr-fe/lorthiois/DeconBenchmark/logs/${DATASET_PREFIX}_${SAMPLE_FILTER}}"

# Create output directory
OUTPUT_BASE_DIR="/work/gr-fe/lorthiois/DeconBenchmark/benchmark_results"
OUTPUT_DIR="${OUTPUT_BASE_DIR}/${DATASET_PREFIX}/benchmarks"
mkdir -p $OUTPUT_DIR

# Ensure log directory exists
mkdir -p "$GLOBAL_LOG_DIR"

module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/
module load r

start=`date +%s`
echo "START AT $(date)"
echo "Processing dataset: $DATASET_PREFIX with filter: $SAMPLE_FILTER"

Rscript ${SCRIPT} ${RLIBRARY} ${DATASET_PREFIX} ${OUTPUT_BASE_DIR} ${INCLUDE_OVERALL_GT} ${SAMPLE_FILTER}

end=`date +%s`
runtime=$((end-start))
echo "Runtime: $runtime seconds"
echo "Benchmark results saved to: $OUTPUT_DIR"
# Record job completion in log file
echo "detailed_metrics_completed:$(date)" >> "$GLOBAL_LOG_DIR/job_status.log"