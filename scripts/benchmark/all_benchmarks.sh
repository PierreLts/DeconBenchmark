#!/bin/bash
# run_benchmarks.sh - Run all benchmarking tasks for a dataset

# Default parameters
DEFAULT_DATASET_PREFIX="PBMC_D100-bulk"
DEFAULT_INCLUDE_OVERALL_GT="TRUE"
DEFAULT_SAMPLE_FILTER="AB"

# Parse command line arguments
DATASET_PREFIX="${1:-$DEFAULT_DATASET_PREFIX}"
INCLUDE_OVERALL_GT="${2:-$DEFAULT_INCLUDE_OVERALL_GT}"
SAMPLE_FILTER="${3:-$DEFAULT_SAMPLE_FILTER}"

# Paths
GLOBAL_LOG_DIR="/work/gr-fe/lorthiois/DeconBenchmark/logs/${DATASET_PREFIX}_${SAMPLE_FILTER}"
mkdir -p "${GLOBAL_LOG_DIR}"
LOG_FILE="${GLOBAL_LOG_DIR}/benchmarking_${DATASET_PREFIX}_${SAMPLE_FILTER}_$(date +%Y%m%d_%H%M%S).log"

SCRIPT_DIR="/work/gr-fe/lorthiois/DeconBenchmark/scripts/benchmark"

echo "===== Starting benchmarking for ${DATASET_PREFIX} =====" | tee -a "${LOG_FILE}"
echo "Include overall ground truth: ${INCLUDE_OVERALL_GT}" | tee -a "${LOG_FILE}"
echo "Sample filter: ${SAMPLE_FILTER}" | tee -a "${LOG_FILE}"
echo "Job mappings will be saved to: ${GLOBAL_LOG_DIR}" | tee -a "${LOG_FILE}"

# 1. Run per-sample plots for all methods
echo "Submitting paired plot generation jobs..." | tee -a "${LOG_FILE}"
PLOT_JOB_ID=$(sbatch "${SCRIPT_DIR}/all_per_sample_plot.sh" "${DATASET_PREFIX}" "${INCLUDE_OVERALL_GT}" "${SAMPLE_FILTER}" "${GLOBAL_LOG_DIR}" | grep -oP 'Submitted batch job \K[0-9]+' || echo "")
echo "Paired plot job ID: ${PLOT_JOB_ID}" | tee -a "${LOG_FILE}"

# 2. Run benchmarking metrics calculation
echo "Submitting benchmarking metrics calculation job..." | tee -a "${LOG_FILE}"
METRICS_JOB_ID=$(sbatch "${SCRIPT_DIR}/detailed_metrics.sh" "${DATASET_PREFIX}" "${INCLUDE_OVERALL_GT}" "${SAMPLE_FILTER}" "${GLOBAL_LOG_DIR}" | grep -oP 'Submitted batch job \K[0-9]+' || echo "")
echo "Metrics job ID: ${METRICS_JOB_ID}" | tee -a "${LOG_FILE}"

# 3. Run benchmarking AB analysis (combined AB samples)
echo "Submitting benchmark summary generation job..." | tee -a "${LOG_FILE}"
SUMMARY_JOB_ID=$(sbatch --dependency=afterok:${METRICS_JOB_ID} "${SCRIPT_DIR}/benchmarking_AB.sh" "/work/gr-fe/R_4.3.1" "${DATASET_PREFIX}" "${SAMPLE_FILTER}" | grep -oP 'Submitted batch job \K[0-9]+' || echo "")
echo "Summary job ID: ${SUMMARY_JOB_ID}" | tee -a "${LOG_FILE}"
echo "summary:${SUMMARY_JOB_ID}" >> "${GLOBAL_LOG_DIR}/stats_job_mapping.txt"

# 4. NEW: Run benchmark_select.sh for specific sample subsets (A, B, and AB)
# These jobs depend on the metrics calculation being completed
for SELECTION in "A" "B" "AB"; do
    echo "Submitting benchmark selection analysis for ${SELECTION} samples..." | tee -a "${LOG_FILE}"
    SELECT_JOB_ID=$(sbatch --dependency=afterok:${METRICS_JOB_ID} "${SCRIPT_DIR}/benchmark_select.sh" "/work/gr-fe/R_4.3.1" "${DATASET_PREFIX}" "${SAMPLE_FILTER}" "${SELECTION}" | grep -oP 'Submitted batch job \K[0-9]+' || echo "")
    echo "Selection (${SELECTION}) job ID: ${SELECT_JOB_ID}" | tee -a "${LOG_FILE}"
    echo "select_${SELECTION}:${SELECT_JOB_ID}" >> "${GLOBAL_LOG_DIR}/stats_job_mapping.txt"
done

# Save job IDs to global stats mapping file
echo "# Stats job mapping for ${DATASET_PREFIX}_${SAMPLE_FILTER} - $(date)" > "${GLOBAL_LOG_DIR}/stats_job_mapping.txt"
echo "paired_plots:${PLOT_JOB_ID}" >> "${GLOBAL_LOG_DIR}/stats_job_mapping.txt"
echo "metrics:${METRICS_JOB_ID}" >> "${GLOBAL_LOG_DIR}/stats_job_mapping.txt"

echo "All benchmarking jobs submitted." | tee -a "${LOG_FILE}"
echo "Check individual job logs for progress." | tee -a "${LOG_FILE}"
echo "Job mapping saved to: ${GLOBAL_LOG_DIR}/stats_job_mapping.txt" | tee -a "${LOG_FILE}"