#!/bin/sh
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 16G
#SBATCH --time 00:30:00
#SBATCH --output=/scratch/lorthiois/logs/%A.o
#SBATCH --error=/scratch/lorthiois/logs/%A.e
#SBATCH --job-name=benchmark_stats

# Default parameters
RLIBRARY="${1:-/work/gr-fe/R_4.3.1}"
DATASET_PREFIX="${2:-TB1}"
SAMPLE_FILTER="${3:-AB}"
# Set benchmark directory to include benchmarks subfolder
BENCHMARK_DIR="${4:-/work/gr-fe/lorthiois/DeconBenchmark/benchmark_results/${DATASET_PREFIX}/benchmarks}"
OUTPUT_DIR="${BENCHMARK_DIR}"  # Use the benchmarks subdirectory directly

# Path to R script
SCRIPT="/work/gr-fe/lorthiois/DeconBenchmark/scripts/benchmark/benchmarking.R"

# Path to detailed metrics CSV - update with filter in filename
METRICS_CSV="${BENCHMARK_DIR}/${DATASET_PREFIX}_detailed_metrics_${SAMPLE_FILTER}.csv"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Load R module
module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/
module load r

start=`date +%s`
echo "START AT $(date)"

# Check if metrics CSV exists
if [ ! -f "$METRICS_CSV" ]; then
    echo "ERROR: Metrics CSV not found at: $METRICS_CSV"
    exit 1
fi

echo "Processing benchmark statistics for $DATASET_PREFIX with filter: $SAMPLE_FILTER"
echo "Using metrics from: $METRICS_CSV"
echo "Output directory: $OUTPUT_DIR"

# Run the R script - pass sample filter parameter
Rscript ${SCRIPT} ${RLIBRARY} ${DATASET_PREFIX} ${METRICS_CSV} ${OUTPUT_DIR} ${SAMPLE_FILTER}

# Check if successful
if [ $? -ne 0 ]; then
    echo "ERROR: R script execution failed"
    exit 1
fi

end=`date +%s`
runtime=$((end-start))
echo "Runtime: $runtime seconds"
echo "Benchmark statistics completed for: $DATASET_PREFIX with filter: $SAMPLE_FILTER"