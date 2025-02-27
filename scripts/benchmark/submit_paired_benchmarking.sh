#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --job-name=submit_paired
#SBATCH --output=/scratch/lorthiois/logs/%A.o
#SBATCH --error=/scratch/lorthiois/logs/%A.e

SCRIPT_DIR="/scratch/lorthiois/scripts"
mkdir -p "$SCRIPT_DIR"
TEMPLATE_DIR="/work/gr-fe/lorthiois/DeconBenchmark/scripts/benchmark"
DECONV_RESULTS_DIR="/work/gr-fe/lorthiois/DeconBenchmark/deconv_results"
BENCHMARK_DIR="/work/gr-fe/lorthiois/DeconBenchmark/paired_benchmark_results"
PER_SAMPLE_GT="/work/gr-fe/lorthiois/DeconBenchmark/generated_data/per_sample_ground_truth_proportions.rda"
GENERAL_GT="/work/gr-fe/lorthiois/DeconBenchmark/generated_data/ground_truth_proportions.rda"
LOG_DIR="/work/gr-fe/lorthiois/DeconBenchmark/logs"

mkdir -p "$BENCHMARK_DIR"
mkdir -p "$LOG_DIR"

# Create mapping files for job tracking
PAIRED_MAPPING_FILE="$LOG_DIR/paired_job_mapping.txt"
PAIRED_PLOT_MAPPING_FILE="$LOG_DIR/paired_plot_job_mapping.txt"

# Clear any existing mapping files
> "$PAIRED_MAPPING_FILE"
> "$PAIRED_PLOT_MAPPING_FILE"

# Log file
SUBMIT_LOG="$LOG_DIR/paired_submission_$(date +%Y%m%d_%H%M%S).log"
echo "==== Starting paired benchmarking submission $(date) ====" | tee -a "$SUBMIT_LOG"

# Find all deconvolution results
RESULTS_FILES=$(find "$DECONV_RESULTS_DIR" -name "results_*.rda")

# Track job IDs for dependency
declare -a JOB_IDS

for RESULTS_FILE in $RESULTS_FILES; do
    # Extract method and dataset name
    FILENAME=$(basename "$RESULTS_FILE")
    METHOD=$(echo "$FILENAME" | sed -E 's/results_([^_]+)_.+\.rda/\1/')
    DATA=$(echo "$FILENAME" | sed -E 's/results_[^_]+_(.+)\.rda/\1/')
    
    echo "Processing $METHOD for dataset $DATA" | tee -a "$SUBMIT_LOG"
    
    # Create method output directory
    METHOD_DIR="$BENCHMARK_DIR/${METHOD}_${DATA}"
    mkdir -p "$METHOD_DIR"
    
    # Create temporary script by copying paired_benchmark.sh and modifying it
    TEMP_SCRIPT="$SCRIPT_DIR/temp_${METHOD}_paired_bench.sh"
    cp "$TEMPLATE_DIR/paired_benchmark.sh" "$TEMP_SCRIPT"
    
    # Modify parameters using sed
    sed -i "s|PER_SAMPLE_GT=.*|PER_SAMPLE_GT=\"$PER_SAMPLE_GT\"|" "$TEMP_SCRIPT"
    sed -i "s|GENERAL_GT=.*|GENERAL_GT=\"$GENERAL_GT\"|" "$TEMP_SCRIPT"
    sed -i "s|RESULTS=.*|RESULTS=\"$RESULTS_FILE\"|" "$TEMP_SCRIPT"
    sed -i "s|OUTPUT_DIR=.*|OUTPUT_DIR=\"$METHOD_DIR\"|" "$TEMP_SCRIPT"
    sed -i "s|#SBATCH --job-name=.*|#SBATCH --job-name=${METHOD}_paired|" "$TEMP_SCRIPT"
    
    # Make script executable
    chmod +x "$TEMP_SCRIPT"
    
    # Submit the job
    JOB_ID=$(sbatch "$TEMP_SCRIPT" | grep -o '[0-9]*')
    if [ -z "$JOB_ID" ]; then
        echo "ERROR: Failed to submit job for $METHOD" | tee -a "$SUBMIT_LOG"
        continue
    fi
    
    echo "Submitted job $JOB_ID for $METHOD paired benchmarking" | tee -a "$SUBMIT_LOG"
    
    # Record the job ID
    JOB_IDS+=("$JOB_ID")
    
    # Add to mapping file
    echo "${METHOD}:${JOB_ID}" >> "$PAIRED_MAPPING_FILE"
done

# Create a dependency string for all jobs
JOB_DEPENDENCY=$(IFS=,; echo "${JOB_IDS[*]}")

# Create and submit comparison job with dependency
if [ ! -z "$JOB_DEPENDENCY" ]; then
    echo "Submitting comparison job dependent on all benchmarking jobs" | tee -a "$SUBMIT_LOG"
    
    # Create a temp directory for the comparison script
    mkdir -p "$BENCHMARK_DIR/comparison"
    
    COMPARE_SCRIPT="$SCRIPT_DIR/temp_compare_paired.sh"
    cp "$TEMPLATE_DIR/compare_paired_methods.sh" "$COMPARE_SCRIPT"
    sed -i "s|BENCHMARK_DIR=.*|BENCHMARK_DIR=\"$BENCHMARK_DIR\"|" "$COMPARE_SCRIPT"
    sed -i "s|OUTPUT_DIR=.*|OUTPUT_DIR=\"$BENCHMARK_DIR/comparison\"|" "$COMPARE_SCRIPT"
    sed -i "s|#SBATCH --job-name=.*|#SBATCH --job-name=compare_paired|" "$COMPARE_SCRIPT"
    
    chmod +x "$COMPARE_SCRIPT"
    
    # Submit with dependency
    COMPARE_JOB_ID=$(sbatch --dependency=afterany:$JOB_DEPENDENCY "$COMPARE_SCRIPT" | grep -o '[0-9]*')
    echo "Submitted comparison job $COMPARE_JOB_ID - will run after all benchmarking jobs complete" | tee -a "$SUBMIT_LOG"
    echo "compare_paired:${COMPARE_JOB_ID}" >> "$PAIRED_MAPPING_FILE"
else
    echo "WARNING: No benchmarking jobs were submitted successfully, skipping comparison job" | tee -a "$SUBMIT_LOG"
fi

echo "==== Submission complete $(date) ====" | tee -a "$SUBMIT_LOG"
echo "Job mapping saved to: $PAIRED_MAPPING_FILE" | tee -a "$SUBMIT_LOG"
echo "All paired benchmarking jobs submitted successfully"