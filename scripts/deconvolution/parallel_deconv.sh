#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=8:00:00
#SBATCH --job-name=parallel_submit
#SBATCH --output=/scratch/lorthiois/logs/%A.o
#SBATCH --error=/scratch/lorthiois/logs/%A.e

# Parse arguments
INPUT_RDA_FILE="$1"
IFS=',' read -r -a MODELS <<< "$2"

# Set paths
SCRIPT_DIR="/scratch/lorthiois/scripts"
mkdir -p "$SCRIPT_DIR"  # Ensure directory exists
OUTPUT_DIR="/work/gr-fe/lorthiois/DeconBenchmark/deconv_results"
BENCHMARK_DIR="/work/gr-fe/lorthiois/DeconBenchmark/benchmark_results"
GROUND_TRUTH="/work/gr-fe/lorthiois/DeconBenchmark/generated_data/ground_truth_proportions.rda"

# Create output directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$BENCHMARK_DIR"

# Extract data name from input path
DATA_NAME=$(basename "$INPUT_RDA_FILE" .rda)

# Log file setup
LOG_DIR="/work/gr-fe/lorthiois/DeconBenchmark/logs"
mkdir -p "$LOG_DIR"
MAIN_LOG="$LOG_DIR/parallel_deconv_$(date +%Y%m%d_%H%M%S).log"

echo "==== Starting parallel deconvolution pipeline $(date) ====" | tee -a "$MAIN_LOG"
echo "Input data: $INPUT_RDA_FILE" | tee -a "$MAIN_LOG"
echo "Models to run: ${MODELS[*]}" | tee -a "$MAIN_LOG"
echo "Total models: ${#MODELS[@]}" | tee -a "$MAIN_LOG"
echo "" | tee -a "$MAIN_LOG"

# Array to store job IDs for each model's deconvolution job
declare -a DECONV_JOB_IDS

# Submit deconvolution jobs for all models in parallel
for MODEL in "${MODELS[@]}"; do
    echo "Submitting deconvolution job for model: $MODEL" | tee -a "$MAIN_LOG"
    
    # Create a temporary job script for this model
    TEMP_SCRIPT="$SCRIPT_DIR/temp_${MODEL}_deconv.sh"
    
    # Copy the template and modify parameters
    cat "$SCRIPT_DIR/deconv_run.sh" | \
        sed "s|input_data=.*|input_data=\"$INPUT_RDA_FILE\"|" | \
        sed "s|output_data=.*|output_data=\"$OUTPUT_DIR\"|" | \
        sed "s|deconv_method=.*|deconv_method=\"$MODEL\"|" | \
        sed "s|#SBATCH --job-name=.*|#SBATCH --job-name=${MODEL}_deconv|" | \
        sed "s|#SBATCH --time=.*|#SBATCH --time=48:00:00|" > "$TEMP_SCRIPT"
        
    # Make executable
    chmod +x "$TEMP_SCRIPT"
    
    # Submit job
    JOB_ID=$(sbatch "$TEMP_SCRIPT" | grep -o '[0-9]*')
    if [ -z "$JOB_ID" ]; then
        echo "ERROR: Failed to submit deconvolution job for $MODEL" | tee -a "$MAIN_LOG"
        continue
    fi
    
    echo "Submitted job $JOB_ID for $MODEL deconvolution" | tee -a "$MAIN_LOG"
    DECONV_JOB_IDS+=("$JOB_ID")
    
    # Store the model-job ID mapping for later use
    echo "${MODEL}:${JOB_ID}" >> "$LOG_DIR/model_job_mapping.txt"
done


# Store job dependencies for evaluation
JOB_DEPENDENCY=$(IFS=,; echo "${DECONV_JOB_IDS[*]}")

# Submit evaluation job with dependency
EVAL_SCRIPT="$SCRIPT_DIR/evaluate_results.sh"
echo "Submitting evaluation job dependent on all deconvolution jobs" | tee -a "$MAIN_LOG"
EVAL_JOB_ID=$(sbatch --dependency=afterany:$JOB_DEPENDENCY "$EVAL_SCRIPT" | grep -o '[0-9]*')
echo "Submitted evaluation job $EVAL_JOB_ID - will run after all deconvolution jobs complete" | tee -a "$MAIN_LOG"

echo "All deconvolution jobs submitted successfully. Monitoring job completion." | tee -a "$MAIN_LOG"

# The script can exit here, and you can use a separate script to monitor
# and process the evaluation after deconvolution jobs complete