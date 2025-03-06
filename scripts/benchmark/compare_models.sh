#!/bin/sh
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 32G
#SBATCH --time 1:00:00
#SBATCH --output=/scratch/lorthiois/logs/%A.o
#SBATCH --error=/scratch/lorthiois/logs/%A.e
#SBATCH --job-name=compare_models

set -e
set -x

DATASET_PREFIX="TB"  # Override with dataset prefix
BENCHMARK_DIR="/work/gr-fe/lorthiois/DeconBenchmark/benchmark_results/${DATASET_PREFIX}"
OUTPUT_DIR="${BENCHMARK_DIR}/comparison"

# Create directory for comparison results
mkdir -p "${OUTPUT_DIR}"

# Run paired benchmark stats to generate comparison metrics
RLIBRARY="/work/gr-fe/R_4.3.1"
BENCHMARK_SCRIPT="/work/gr-fe/lorthiois/DeconBenchmark/scripts/benchmark/paired_benchmark_stats.sh"

# Run the benchmark script
bash ${BENCHMARK_SCRIPT} ${DATASET_PREFIX} TRUE

echo "Comparison completed. Results saved to: ${OUTPUT_DIR}"