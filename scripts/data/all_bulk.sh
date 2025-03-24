#!/bin/bash
# all_bulk.sh - Run all bulk data generation tasks

# Default parameters
DEFAULT_PREFIX="TB_D33348"
DEFAULT_RLIBRARY="/work/gr-fe/R_4.3.1"
DEFAULT_SEURAT_FILE="/work/gr-fe/lorthiois/DeconBenchmark/data/merged_batches.rds"
DEFAULT_BULK_FILE="/work/gr-fe/lorthiois/DeconBenchmark/data/cleaned_feature_counts_matrix.csv"
DEFAULT_MAPPING_FILE="/work/gr-fe/lorthiois/DeconBenchmark/data/mart_export.txt"
DEFAULT_OUTPUT_DIR="/work/gr-fe/lorthiois/DeconBenchmark/generated_data"

# Parse command line arguments
PREFIX="${1:-$DEFAULT_PREFIX}"
RLIBRARY="${2:-$DEFAULT_RLIBRARY}"
SEURAT_FILE="${3:-$DEFAULT_SEURAT_FILE}"
BULK_FILE="${4:-$DEFAULT_BULK_FILE}"
MAPPING_FILE="${5:-$DEFAULT_MAPPING_FILE}"
OUTPUT_DIR="${6:-$DEFAULT_OUTPUT_DIR}"

# Script directories
SCRIPT_DIR="/work/gr-fe/lorthiois/DeconBenchmark/scripts/data"
SCRATCH_DIR="/scratch/lorthiois/scripts"
GLOBAL_LOG_DIR="/work/gr-fe/lorthiois/DeconBenchmark/logs/${PREFIX}"

# Create directories
mkdir -p "$SCRATCH_DIR"
mkdir -p "$GLOBAL_LOG_DIR"
mkdir -p "$OUTPUT_DIR/$PREFIX"

# Set up log file 
LOG_FILE="$GLOBAL_LOG_DIR/all_bulk_$(date +%Y%m%d_%H%M%S).log"
JOB_MAPPING_FILE="$GLOBAL_LOG_DIR/bulk_job_mapping.txt"
echo "# Bulk job mapping file for ${PREFIX} - Created $(date)" > "$JOB_MAPPING_FILE"

echo "==== Starting bulk data generation tasks $(date) ====" | tee -a "$LOG_FILE"
echo "PREFIX: $PREFIX" | tee -a "$LOG_FILE"
echo "Output directory: $OUTPUT_DIR/$PREFIX" | tee -a "$LOG_FILE"

# Function to submit jobs
submit_job() {
    local script_name=$1
    local job_name=$2
    local template_path="$SCRIPT_DIR/$script_name"
    local temp_script="$SCRATCH_DIR/${job_name}_${PREFIX}.sh"
    
    # Copy template and modify parameters
    cat "$template_path" | \
        sed "s|RLIBRARY=.*|RLIBRARY=\"$RLIBRARY\"|" | \
        sed "s|INPUT_DATA=.*|INPUT_DATA=\"$BULK_FILE\"|" | \
        sed "s|OUTPUT_BASE_DIR=.*|OUTPUT_BASE_DIR=\"$OUTPUT_DIR\"|" | \
        sed "s|PREFIX=.*|PREFIX=\"$PREFIX\"|" | \
        sed "s|SEURAT_FILE=.*|SEURAT_FILE=\"$SEURAT_FILE\"|" | \
        sed "s|MAPPING_FILE=.*|MAPPING_FILE=\"$MAPPING_FILE\"|" | \
        sed "s|#SBATCH --job-name=.*|#SBATCH --job-name=${job_name}_${PREFIX}|" > "$temp_script"
    
    chmod +x "$temp_script"
    
    # Submit job
    echo "Submitting $job_name job..." | tee -a "$LOG_FILE"
    JOB_ID=$(sbatch "$temp_script" | grep -o '[0-9]*')
    
    if [ -z "$JOB_ID" ]; then
        echo "ERROR: Failed to submit job for $job_name" | tee -a "$LOG_FILE"
        return 1
    fi
    
    echo "${job_name}:${JOB_ID}" >> "$JOB_MAPPING_FILE"
    echo "Submitted job $JOB_ID for $job_name" | tee -a "$LOG_FILE"
    return 0
}

# Submit bulk generation job
submit_job "bulk_generation.sh" "bulk_gen"

# Submit pseudobulk generation job
submit_job "pseudobulk.sh" "pseudobulk_gen"

# Submit bulk randomizer job
submit_job "bulk_randomizer.sh" "bulk_random_gen"

echo "All bulk data generation jobs submitted successfully." | tee -a "$LOG_FILE"
echo "Job mapping saved to: $JOB_MAPPING_FILE" | tee -a "$LOG_FILE"
echo "Check individual job logs for progress." | tee -a "$LOG_FILE"