#!/bin/bash
# run_multiple_models_fixed.sh
# This script runs multiple deconvolution models sequentially,
# followed by evaluation and visualization for each model.
#
# Usage: ./run_multiple_models_fixed.sh INPUT_RDA_FILE MODEL1,MODEL2,...

set -o pipefail  # Exit on pipe failures

# Validate arguments
if [ $# -lt 2 ]; then
  echo "Usage: $0 <input_rda_file> <comma_separated_model_list>"
  echo "Example: $0 /path/to/Batch1.rda MuSiC,CIBERSORT,RNA-Sieve"
  exit 1
fi

# Parse arguments
INPUT_RDA_FILE="$1"
IFS=',' read -r -a MODELS <<< "$2"

# Set paths
SCRIPT_DIR="/work/gr-fe/lorthiois/DeconBenchmark/scripts"
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
MAIN_LOG="$LOG_DIR/run_multiple_models_$(date +%Y%m%d_%H%M%S).log"

echo "==== Starting deconvolution pipeline $(date) ====" | tee -a "$MAIN_LOG"
echo "Input data: $INPUT_RDA_FILE" | tee -a "$MAIN_LOG"
echo "Models to run: ${MODELS[*]}" | tee -a "$MAIN_LOG"
echo "Total models: ${#MODELS[@]}" | tee -a "$MAIN_LOG"
echo "" | tee -a "$MAIN_LOG"

# Run each model sequentially
for ((i=0; i<${#MODELS[@]}; i++)); do
  MODEL="${MODELS[$i]}"
  MODEL_LOG="$LOG_DIR/model_${MODEL}_$(date +%Y%m%d_%H%M%S).log"
  
  echo "[$((i+1))/${#MODELS[@]}] Processing model: $MODEL ($(date))" | tee -a "$MAIN_LOG"
  
  # 1. Run deconvolution
  echo "  Running deconvolution..." | tee -a "$MAIN_LOG"
  RESULTS_FILE="$OUTPUT_DIR/results_${MODEL}_${DATA_NAME}.rda"
  
  # Modify deconv_run.sh parameters
  cat "$SCRIPT_DIR/deconv_run.sh" | \
    sed "s|input_data=.*|input_data=\"$INPUT_RDA_FILE\"|" | \
    sed "s|output_data=.*|output_data=\"$OUTPUT_DIR\"|" | \
    sed "s|deconv_method=.*|deconv_method=\"$MODEL\"|" | \
    sed "s|#SBATCH --job-name=.*|#SBATCH --job-name=${MODEL}_deconv|" | \
    sed "s|#SBATCH --time=.*|#SBATCH --time=24:00:00|" > \
    "$SCRIPT_DIR/temp_deconv_run.sh"
  
  # Ensure time limit is explicitly set by adding it if not present
  if ! grep -q "#SBATCH --time=" "$SCRIPT_DIR/temp_deconv_run.sh"; then
    sed -i '/#SBATCH --mem/a #SBATCH --time=24:00:00' "$SCRIPT_DIR/temp_deconv_run.sh"
  fi
  
  # Make executable
  chmod +x "$SCRIPT_DIR/temp_deconv_run.sh"
  
  # Submit job and wait for completion
  JOB_ID=$(sbatch "$SCRIPT_DIR/temp_deconv_run.sh" | grep -o '[0-9]*')
  if [ -z "$JOB_ID" ]; then
    echo "  ERROR: Failed to submit deconvolution job for $MODEL" | tee -a "$MAIN_LOG" "$MODEL_LOG"
    continue  # Skip to next model
  fi
  
  echo "  Submitted job $JOB_ID for deconvolution" | tee -a "$MAIN_LOG" "$MODEL_LOG"
  
  # Wait for job completion with timeout check
  MAX_WAIT_TIME=86400  # 24 hours in seconds
  START_TIME=$(date +%s)
  
  while squeue -j $JOB_ID &>/dev/null; do
    CURRENT_TIME=$(date +%s)
    ELAPSED_TIME=$((CURRENT_TIME - START_TIME))
    
    # If job has been running for too long, cancel and move on
    if [ $ELAPSED_TIME -gt $MAX_WAIT_TIME ]; then
      echo "  WARNING: Job $JOB_ID for $MODEL has been running for more than 24 hours. Cancelling." | tee -a "$MAIN_LOG" "$MODEL_LOG"
      scancel $JOB_ID
      break
    fi
    
    sleep 120  # Check every 2 minutes
  done
  
  # Check if results file exists
  if [ ! -f "$RESULTS_FILE" ]; then
    echo "  ERROR: Deconvolution failed for $MODEL - results file not created" | tee -a "$MAIN_LOG" "$MODEL_LOG"
    continue  # Skip to next model
  fi
  
  # 2. Run model statistics
  echo "  Running benchmarking..." | tee -a "$MAIN_LOG"
  
  # Modify model_stats.sh parameters
  cat "$SCRIPT_DIR/model_stats.sh" | \
    sed "s|GROUND_TRUTH=.*|GROUND_TRUTH=\"$GROUND_TRUTH\"|" | \
    sed "s|RESULTS=.*|RESULTS=\"$RESULTS_FILE\"|" | \
    sed "s|OUTPUT_DIR=.*|OUTPUT_DIR=\"$BENCHMARK_DIR\"|" | \
    sed "s|#SBATCH --job-name=.*|#SBATCH --job-name=${MODEL}_stats|" | \
    sed "s|#SBATCH --time=.*|#SBATCH --time=2:00:00|" > \
    "$SCRIPT_DIR/temp_model_stats.sh"
  
  # Ensure time limit is explicitly set
  if ! grep -q "#SBATCH --time=" "$SCRIPT_DIR/temp_model_stats.sh"; then
    sed -i '/#SBATCH --mem/a #SBATCH --time=2:00:00' "$SCRIPT_DIR/temp_model_stats.sh"
  fi
  
  # Make executable
  chmod +x "$SCRIPT_DIR/temp_model_stats.sh"
  
  # Submit job and wait for completion
  JOB_ID=$(sbatch "$SCRIPT_DIR/temp_model_stats.sh" | grep -o '[0-9]*')
  if [ -z "$JOB_ID" ]; then
    echo "  ERROR: Failed to submit benchmarking job for $MODEL" | tee -a "$MAIN_LOG" "$MODEL_LOG"
    continue  # Skip visualization but proceed to next model
  fi
  
  echo "  Submitted job $JOB_ID for benchmarking" | tee -a "$MAIN_LOG" "$MODEL_LOG"
  
  # Wait for job completion
  while squeue -j $JOB_ID &>/dev/null; do
    sleep 30
  done
  
  # 3. Run visualization
  echo "  Running visualization..." | tee -a "$MAIN_LOG"
  
  # Modify cell_proportion_plot.sh parameters
  cat "$SCRIPT_DIR/cell_proportion_plot.sh" | \
    sed "s|GROUND_TRUTH=.*|GROUND_TRUTH=\"$GROUND_TRUTH\"|" | \
    sed "s|RESULTS=.*|RESULTS=\"$RESULTS_FILE\"|" | \
    sed "s|OUTPUT_DIR=.*|OUTPUT_DIR=\"$BENCHMARK_DIR\"|" | \
    sed "s|#SBATCH --job-name=.*|#SBATCH --job-name=${MODEL}_plot|" | \
    sed "s|#SBATCH --time=.*|#SBATCH --time=1:00:00|" > \
    "$SCRIPT_DIR/temp_cell_proportion_plot.sh"
  
  # Ensure time limit is explicitly set
  if ! grep -q "#SBATCH --time=" "$SCRIPT_DIR/temp_cell_proportion_plot.sh"; then
    sed -i '/#SBATCH --mem/a #SBATCH --time=1:00:00' "$SCRIPT_DIR/temp_cell_proportion_plot.sh"
  fi
  
  # Make executable
  chmod +x "$SCRIPT_DIR/temp_cell_proportion_plot.sh"
  
  # Submit job and wait for completion
  JOB_ID=$(sbatch "$SCRIPT_DIR/temp_cell_proportion_plot.sh" | grep -o '[0-9]*')
  if [ -z "$JOB_ID" ]; then
    echo "  ERROR: Failed to submit visualization job for $MODEL" | tee -a "$MAIN_LOG" "$MODEL_LOG"
    continue  # Skip to next model
  fi
  
  echo "  Submitted job $JOB_ID for visualization" | tee -a "$MAIN_LOG" "$MODEL_LOG"
  
  # Wait for job completion
  while squeue -j $JOB_ID &>/dev/null; do
    sleep 30
  done
  
  echo "  Completed processing for $MODEL" | tee -a "$MAIN_LOG"
done

# Run comparison script after all models are processed
if [ ${#MODELS[@]} -gt 1 ]; then
  echo "Running model comparison..." | tee -a "$MAIN_LOG"
  
  # Modify compare_models.sh parameters
  cat "$SCRIPT_DIR/compare_models.sh" | \
    sed "s|BENCHMARK_DIR=.*|BENCHMARK_DIR=\"$BENCHMARK_DIR\"|" | \
    sed "s|OUTPUT_DIR=.*|OUTPUT_DIR=\"$BENCHMARK_DIR/comparison\"|" | \
    sed "s|#SBATCH --job-name=.*|#SBATCH --job-name=compare_models|" | \
    sed "s|#SBATCH --time=.*|#SBATCH --time=1:00:00|" > \
    "$SCRIPT_DIR/temp_compare_models.sh"
  
  # Ensure time limit is explicitly set
  if ! grep -q "#SBATCH --time=" "$SCRIPT_DIR/temp_compare_models.sh"; then
    sed -i '/#SBATCH --mem/a #SBATCH --time=1:00:00' "$SCRIPT_DIR/temp_compare_models.sh"
  fi
  
  # Make executable
  chmod +x "$SCRIPT_DIR/temp_compare_models.sh"
  
  # Submit job and wait for completion
  JOB_ID=$(sbatch "$SCRIPT_DIR/temp_compare_models.sh" | grep -o '[0-9]*')
  if [ -z "$JOB_ID" ]; then
    echo "ERROR: Failed to submit model comparison job" | tee -a "$MAIN_LOG"
  else
    echo "Submitted job $JOB_ID for model comparison" | tee -a "$MAIN_LOG"
    
    # Wait for job completion
    while squeue -j $JOB_ID &>/dev/null; do
      sleep 30
    done
  fi
fi

# Cleanup temporary files
rm -f "$SCRIPT_DIR/temp_deconv_run.sh" "$SCRIPT_DIR/temp_model_stats.sh" \
      "$SCRIPT_DIR/temp_cell_proportion_plot.sh" "$SCRIPT_DIR/temp_compare_models.sh"

echo "==== Finished deconvolution pipeline $(date) ====" | tee -a "$MAIN_LOG"
echo "See detailed logs in $LOG_DIR" | tee -a "$MAIN_LOG"