#!/bin/bash
# run_benchmarks.sh - Run all benchmarking tasks for a dataset
#
# Usage: ./run_benchmarks.sh [DATASET_PREFIX] [INCLUDE_OVERALL_GT]

# Default parameters
DEFAULT_DATASET_PREFIX="TB1"
DEFAULT_INCLUDE_OVERALL_GT="TRUE"

# Parse command line arguments
DATASET_PREFIX="${1:-$DEFAULT_DATASET_PREFIX}"
INCLUDE_OVERALL_GT="${2:-$DEFAULT_INCLUDE_OVERALL_GT}"

# Paths
SCRIPT_DIR="/work/gr-fe/lorthiois/DeconBenchmark/scripts/benchmark"
LOG_DIR="/work/gr-fe/lorthiois/DeconBenchmark/logs/${DATASET_PREFIX}"
mkdir -p "${LOG_DIR}"

# Log file
LOG_FILE="${LOG_DIR}/benchmarking_${DATASET_PREFIX}_$(date +%Y%m%d_%H%M%S).log"

echo "===== Starting benchmarking for ${DATASET_PREFIX} =====" | tee -a "${LOG_FILE}"
echo "Include overall ground truth: ${INCLUDE_OVERALL_GT}" | tee -a "${LOG_FILE}"

# 1. Run per-sample plots for all methods
echo "Submitting paired plot generation jobs..." | tee -a "${LOG_FILE}"
sbatch "${SCRIPT_DIR}/per_sample_multi_plot.sh" "${DATASET_PREFIX}" "${INCLUDE_OVERALL_GT}" | tee -a "${LOG_FILE}"

# 2. Run benchmarking metrics calculation
echo "Submitting benchmarking metrics calculation job..." | tee -a "${LOG_FILE}"
sbatch "${SCRIPT_DIR}/paired_benchmark_stats.sh" "${DATASET_PREFIX}" "${INCLUDE_OVERALL_GT}" | tee -a "${LOG_FILE}"

echo "All benchmarking jobs submitted." | tee -a "${LOG_FILE}"
echo "Check individual job logs for progress."