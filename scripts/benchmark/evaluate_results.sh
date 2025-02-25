#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --job-name=eval_submit
#SBATCH --output=/scratch/lorthiois/logs/%A.o
#SBATCH --error=/scratch/lorthiois/logs/%A.e

# evaluate_results.sh - Run evaluation for completed deconvolution jobs
# 
# Usage: ./evaluate_results.sh

SCRIPT_DIR="/scratch/lorthiois/scripts"
mkdir -p "$SCRIPT_DIR"
TEMPLATE_DIR="/work/gr-fe/lorthiois/DeconBenchmark/scripts/benchmark"
OUTPUT_DIR="/work/gr-fe/lorthiois/DeconBenchmark/deconv_results"
BENCHMARK_DIR="/work/gr-fe/lorthiois/DeconBenchmark/benchmark_results"
GROUND_TRUTH="/work/gr-fe/lorthiois/DeconBenchmark/generated_data/ground_truth_proportions.rda"
LOG_DIR="/work/gr-fe/lorthiois/DeconBenchmark/logs"
MAPPING_FILE="$LOG_DIR/model_job_mapping.txt"

DATA_NAME=$(basename $(ls $OUTPUT_DIR/results_*_*.rda | head -1) | sed -E 's/results_[^_]+_(.+)\.rda/\1/')

# Log file
EVAL_LOG="$LOG_DIR/evaluation_$(date +%Y%m%d_%H%M%S).log"

echo "==== Starting evaluation pipeline $(date) ====" | tee -a "$EVAL_LOG"

# Read the model-job mapping file
if [ ! -f "$MAPPING_FILE" ]; then
    echo "ERROR: Mapping file not found: $MAPPING_FILE" | tee -a "$EVAL_LOG"
    exit 1
fi

# For each model, check if result file exists and run evaluation
while IFS=: read -r MODEL JOB_ID; do
    RESULTS_FILE="$OUTPUT_DIR/results_${MODEL}_${DATA_NAME}.rda"
    
    # Check if results file exists
    if [ ! -f "$RESULTS_FILE" ]; then
        echo "Results not found for $MODEL - skipping evaluation" | tee -a "$EVAL_LOG"
        continue
    fi
    
    echo "Processing evaluation for $MODEL" | tee -a "$EVAL_LOG"
    
    # Run model statistics
    STATS_SCRIPT="$SCRIPT_DIR/temp_${MODEL}_stats.sh"
    cat "$TEMPLATE_DIR/model_stats.sh" | \
        sed "s|GROUND_TRUTH=.*|GROUND_TRUTH=\"$GROUND_TRUTH\"|" | \
        sed "s|RESULTS=.*|RESULTS=\"$RESULTS_FILE\"|" | \
        sed "s|OUTPUT_DIR=.*|OUTPUT_DIR=\"$BENCHMARK_DIR\"|" | \
        sed "s|#SBATCH --job-name=.*|#SBATCH --job-name=${MODEL}_stats|" | \
        sed "s|#SBATCH --time=.*|#SBATCH --time=2:00:00|" > "$STATS_SCRIPT"
    
    chmod +x "$STATS_SCRIPT"
    
    # Submit statistics job
    STATS_JOB_ID=$(sbatch "$STATS_SCRIPT" | grep -o '[0-9]*')
    echo "Submitted job $STATS_JOB_ID for $MODEL statistics" | tee -a "$EVAL_LOG"
    
    # Run visualization (with dependency on stats job)
    PLOT_SCRIPT="$SCRIPT_DIR/temp_${MODEL}_plot.sh"
    cat "$SCRIPT_DIR/cell_proportion_plot.sh" | \
        sed "s|GROUND_TRUTH=.*|GROUND_TRUTH=\"$GROUND_TRUTH\"|" | \
        sed "s|RESULTS=.*|RESULTS=\"$RESULTS_FILE\"|" | \
        sed "s|OUTPUT_DIR=.*|OUTPUT_DIR=\"$BENCHMARK_DIR\"|" | \
        sed "s|#SBATCH --job-name=.*|#SBATCH --job-name=${MODEL}_plot|" | \
        sed "s|#SBATCH --time=.*|#SBATCH --time=1:00:00|" > "$PLOT_SCRIPT"
    
    chmod +x "$PLOT_SCRIPT"
    
    # Submit plot job with dependency on stats job
    PLOT_JOB_ID=$(sbatch --dependency=afterok:$STATS_JOB_ID "$PLOT_SCRIPT" | grep -o '[0-9]*')
    echo "Submitted job $PLOT_JOB_ID for $MODEL visualization" | tee -a "$EVAL_LOG"
    
done < "$MAPPING_FILE"

# Submit comparison job
echo "Running model comparison..." | tee -a "$EVAL_LOG"
COMPARE_SCRIPT="$SCRIPT_DIR/temp_compare_models.sh"
cat "$SCRIPT_DIR/compare_models.sh" | \
    sed "s|BENCHMARK_DIR=.*|BENCHMARK_DIR=\"$BENCHMARK_DIR\"|" | \
    sed "s|OUTPUT_DIR=.*|OUTPUT_DIR=\"$BENCHMARK_DIR/comparison\"|" | \
    sed "s|#SBATCH --job-name=.*|#SBATCH --job-name=compare_models|" | \
    sed "s|#SBATCH --time=.*|#SBATCH --time=2:00:00|" > "$COMPARE_SCRIPT"

chmod +x "$COMPARE_SCRIPT"
COMPARE_JOB_ID=$(sbatch "$COMPARE_SCRIPT" | grep -o '[0-9]*')
echo "Submitted job $COMPARE_JOB_ID for model comparison" | tee -a "$EVAL_LOG"

echo "==== Evaluation pipeline completed $(date) ====" | tee -a "$EVAL_LOG"