#!/bin/bash
# run_benchmarks.sh - Run all benchmarking tasks for a dataset
#
# Usage: ./run_benchmarks.sh [DATASET_PREFIX] [INCLUDE_OVERALL_GT] [SAMPLE_FILTER]

# Default parameters
DEFAULT_DATASET_PREFIX="TB1"
DEFAULT_INCLUDE_OVERALL_GT="TRUE"
DEFAULT_SAMPLE_FILTER="A"

# Parse command line arguments
DATASET_PREFIX="${1:-$DEFAULT_DATASET_PREFIX}"
INCLUDE_OVERALL_GT="${2:-$DEFAULT_INCLUDE_OVERALL_GT}"
SAMPLE_FILTER="${3:-$DEFAULT_SAMPLE_FILTER}"

# Paths
GLOBAL_LOG_DIR="/work/gr-fe/lorthiois/DeconBenchmark/logs/${DATASET_PREFIX}"
mkdir -p "${GLOBAL_LOG_DIR}"
LOG_FILE="${GLOBAL_LOG_DIR}/benchmarking_${DATASET_PREFIX}_${SAMPLE_FILTER}_$(date +%Y%m%d_%H%M%S).log"

SCRIPT_DIR="/work/gr-fe/lorthiois/DeconBenchmark/scripts/benchmark"
LOG_DIR="/work/gr-fe/lorthiois/DeconBenchmark/logs/${DATASET_PREFIX}"
mkdir -p "${LOG_DIR}"

echo "===== Starting benchmarking for ${DATASET_PREFIX} =====" | tee -a "${LOG_FILE}"
echo "Include overall ground truth: ${INCLUDE_OVERALL_GT}" | tee -a "${LOG_FILE}"
echo "Sample filter: ${SAMPLE_FILTER}" | tee -a "${LOG_FILE}"
echo "Job mappings will be saved to: ${GLOBAL_LOG_DIR}" | tee -a "${LOG_FILE}"

# 1. Run per-sample plots for all methods
echo "Submitting paired plot generation jobs..." | tee -a "${LOG_FILE}"
sbatch "${SCRIPT_DIR}/per_sample_multi_plot.sh" "${DATASET_PREFIX}" "${INCLUDE_OVERALL_GT}" "${SAMPLE_FILTER}" "${GLOBAL_LOG_DIR}" | tee -a "${LOG_FILE}"
# 2. Run benchmarking metrics calculation
echo "Submitting benchmarking metrics calculation job..." | tee -a "${LOG_FILE}"
sbatch "${SCRIPT_DIR}/detailed_metrics.sh" "${DATASET_PREFIX}" "${INCLUDE_OVERALL_GT}" "${SAMPLE_FILTER}" "${GLOBAL_LOG_DIR}" | tee -a "${LOG_FILE}"

echo "All benchmarking jobs submitted." | tee -a "${LOG_FILE}"
echo "Check individual job logs for progress."