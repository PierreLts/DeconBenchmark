#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --job-name=paired_plots
#SBATCH --output=/scratch/lorthiois/logs/%A.o
#SBATCH --error=/scratch/lorthiois/logs/%A.e

# Default parameters
DATASET_PREFIX="${1:-TB}"
INCLUDE_OVERALL_GT="${2:-TRUE}"
SAMPLE_FILTER="${3:-AB}"  # New parameter for sample filter
GLOBAL_LOG_DIR="${4:-/work/gr-fe/lorthiois/DeconBenchmark/logs/${DATASET_PREFIX}_${SAMPLE_FILTER}}"

# Set paths
SCRIPT_DIR="/scratch/lorthiois/scripts"
mkdir -p "$SCRIPT_DIR"
TEMPLATE_DIR="/work/gr-fe/lorthiois/DeconBenchmark/scripts/benchmark"
DECONV_RESULTS_DIR="/work/gr-fe/lorthiois/DeconBenchmark/deconv_results/${DATASET_PREFIX}"
BENCHMARK_DIR="/work/gr-fe/lorthiois/DeconBenchmark/benchmark_results/${DATASET_PREFIX}"

# Create directories
mkdir -p "${BENCHMARK_DIR}"
mkdir -p "$GLOBAL_LOG_DIR"

# Start with clean mapping file
MAPPING_FILE="$GLOBAL_LOG_DIR/paired_plot_job_mapping.txt"
> "$MAPPING_FILE"

# Find all result files for this dataset with the correct filter
RESULT_FILES=$(find $DECONV_RESULTS_DIR -name "results_*_${SAMPLE_FILTER}.rda")

# Check if we found any files
if [ -z "$RESULT_FILES" ]; then
    echo "No result files found in $DECONV_RESULTS_DIR matching filter: $SAMPLE_FILTER"
    exit 1
fi

echo "Found $(echo "$RESULT_FILES" | wc -l) result files to process with filter: $SAMPLE_FILTER"

# Process each result file
for RESULTS_FILE in $RESULT_FILES; do
    # Extract method name from filename, accounting for filter suffix
    METHOD=$(basename $RESULTS_FILE | sed -E "s/results_([^_]+)_${SAMPLE_FILTER}\.rda/\1/")
    
    echo "Processing $METHOD..."
    
    # Create temporary script for this method
    TEMP_SCRIPT="$SCRIPT_DIR/temp_${METHOD}_${DATASET_PREFIX}_${SAMPLE_FILTER}_paired_plot.sh"
    
    # Copy the template and modify parameters
    cat "$TEMPLATE_DIR/per_sample_cell_proportion_plot.sh" | \
        sed "s|DATASET_PREFIX=.*|DATASET_PREFIX=\"$DATASET_PREFIX\"|" | \
        sed "s|METHOD=.*|METHOD=\"$METHOD\"|" | \
        sed "s|SAMPLE_FILTER=.*|SAMPLE_FILTER=\"$SAMPLE_FILTER\"|" | \
        sed "s|INCLUDE_OVERALL_GT=.*|INCLUDE_OVERALL_GT=\"$INCLUDE_OVERALL_GT\"|" | \
        sed "s|#SBATCH --job-name=.*|#SBATCH --job-name=${METHOD}_${DATASET_PREFIX}_${SAMPLE_FILTER}_plot|" > "$TEMP_SCRIPT"
        
    chmod +x "$TEMP_SCRIPT"
    
    # Submit job
    JOB_ID=$(sbatch "$TEMP_SCRIPT" | grep -o '[0-9]*')
    
    if [ -z "$JOB_ID" ]; then
        echo "ERROR: Failed to submit job for $METHOD"
        continue
    fi
    
    echo "Submitted job $JOB_ID for $METHOD paired plot"
    # Record mapping
    echo "${METHOD}:${JOB_ID}" >> "$MAPPING_FILE"
done

echo "All paired plot jobs submitted. Mapping saved to $MAPPING_FILE"

# Submit a comparison job to run after all plots are created
COMPARE_SCRIPT="$SCRIPT_DIR/temp_compare_${DATASET_PREFIX}_${SAMPLE_FILTER}_paired_plots.sh"
cat "$TEMPLATE_DIR/compare_models.sh" | \
    sed "s|DATASET_PREFIX=.*|DATASET_PREFIX=\"$DATASET_PREFIX\"|" | \
    sed "s|BENCHMARK_DIR=.*|BENCHMARK_DIR=\"$BENCHMARK_DIR\"|" | \
    sed "s|OUTPUT_DIR=.*|OUTPUT_DIR=\"$BENCHMARK_DIR/comparison_${SAMPLE_FILTER}\"|" | \
    sed "s|#SBATCH --job-name=.*|#SBATCH --job-name=compare_${DATASET_PREFIX}_${SAMPLE_FILTER}|" > "$COMPARE_SCRIPT"

chmod +x "$COMPARE_SCRIPT"

# Get all job IDs for dependency
JOB_IDS=$(cut -d':' -f2 "$MAPPING_FILE" | paste -sd "," -)

if [ -n "$JOB_IDS" ]; then
    COMPARE_JOB_ID=$(sbatch --dependency=afterany:$JOB_IDS "$COMPARE_SCRIPT" | grep -o '[0-9]*')
    echo "Submitted comparison job $COMPARE_JOB_ID, dependent on all plot jobs"
    
    # Add comparison job to mapping file
    echo "compare_${DATASET_PREFIX}_${SAMPLE_FILTER}:${COMPARE_JOB_ID}" >> "$MAPPING_FILE"
else
    echo "No plot jobs were submitted successfully, skipping comparison job"
fi