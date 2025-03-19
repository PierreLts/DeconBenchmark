#!/bin/bash
# master_data_generation.sh - Run all data generation tasks in parallel

# Default parameters
DEFAULT_SEURAT_FILE="/work/gr-fe/lorthiois/DeconBenchmark/data/GFB-33245_HFKJMDSXC_2_scRNAseqWTATBseverityrun1_Seurat.rds"
DEFAULT_BULK_FILE="/work/gr-fe/lorthiois/DeconBenchmark/data/cleaned_feature_counts_matrix.csv"
DEFAULT_MAPPING_FILE="/work/gr-fe/lorthiois/DeconBenchmark/data/mart_export.txt"
DEFAULT_OUTPUT_DIR="/work/gr-fe/lorthiois/DeconBenchmark/generated_data"
DEFAULT_PREFIX="TB1"
DEFAULT_RLIBRARY="/work/gr-fe/R_4.3.1"

# Parse command line arguments
SEURAT_FILE=${1:-$DEFAULT_SEURAT_FILE}
BULK_FILE=${2:-$DEFAULT_BULK_FILE}
MAPPING_FILE=${3:-$DEFAULT_MAPPING_FILE}
OUTPUT_DIR=${4:-$DEFAULT_OUTPUT_DIR}
PREFIX=${5:-$DEFAULT_PREFIX}
RLIBRARY=${6:-$DEFAULT_RLIBRARY}

# Create prefix-specific subdirectory for output data
SUBDIR_PATH="$OUTPUT_DIR/$PREFIX"
mkdir -p $SUBDIR_PATH

# Set up scratch log directory
SCRATCH_LOG_DIR="/scratch/lorthiois/logs"
mkdir -p $SCRATCH_LOG_DIR

# Create prefix-specific subdirectory for temporary scripts
SCRATCH_SCRIPT_DIR="/scratch/lorthiois/scripts/${PREFIX}"
mkdir -p $SCRATCH_SCRIPT_DIR

# Set up log file in output directory for logging the script execution
OUTPUT_LOG_DIR="$SUBDIR_PATH/logs"
mkdir -p $OUTPUT_LOG_DIR
MAIN_LOG="$OUTPUT_LOG_DIR/master_data_generation_$(date +%Y%m%d_%H%M%S).log"

echo "==== Starting master data generation $(date) ====" | tee -a "$MAIN_LOG"
echo "Seurat file: $SEURAT_FILE" | tee -a "$MAIN_LOG"
echo "Bulk file: $BULK_FILE" | tee -a "$MAIN_LOG" 
echo "Mapping file: $MAPPING_FILE" | tee -a "$MAIN_LOG"
echo "Output directory: $SUBDIR_PATH" | tee -a "$MAIN_LOG"
echo "File prefix: $PREFIX" | tee -a "$MAIN_LOG"
echo "R library path: $RLIBRARY" | tee -a "$MAIN_LOG"
echo "Scratch log directory: $SCRATCH_LOG_DIR" | tee -a "$MAIN_LOG"
echo "Scratch script directory: $SCRATCH_SCRIPT_DIR" | tee -a "$MAIN_LOG"
echo "" | tee -a "$MAIN_LOG"

# Base directory where scripts are located
SCRIPT_DIR="/work/gr-fe/lorthiois/DeconBenchmark/scripts/data"

# Function to create temp script and submit job
submit_job() {
    local job_name=$1
    local script_content=$2
    local cpus=$3
    local mem=$4
    local time=$5
    local dependency=$6
    
    # Create temporary script
    local temp_script="$SCRATCH_SCRIPT_DIR/${job_name}_${PREFIX}.sh"
    echo "#!/bin/bash" > "$temp_script"
    echo "#SBATCH --job-name=${job_name}_${PREFIX}" >> "$temp_script"
    echo "#SBATCH --output=${SCRATCH_LOG_DIR}/%j.o" >> "$temp_script"
    echo "#SBATCH --error=${SCRATCH_LOG_DIR}/%j.e" >> "$temp_script"
    echo "#SBATCH --nodes=1" >> "$temp_script"
    echo "#SBATCH --ntasks=1" >> "$temp_script"
    echo "#SBATCH --cpus-per-task=${cpus}" >> "$temp_script"
    echo "#SBATCH --mem=${mem}" >> "$temp_script"
    echo "#SBATCH --time=${time}" >> "$temp_script"
    
    # Add dependency if provided
    if [ ! -z "$dependency" ]; then
        echo "#SBATCH --dependency=afterok:${dependency}" >> "$temp_script"
    fi
    
    # Add module loading and script execution
    echo "module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/" >> "$temp_script"
    echo "module load r" >> "$temp_script"
    echo "" >> "$temp_script"
    echo "$script_content" >> "$temp_script"
    
    # Make script executable
    chmod +x "$temp_script"
    
    # Submit job and capture output with proper error handling
    local sbatch_output=$(sbatch "$temp_script" 2>&1)
    local sbatch_status=$?
    
    if [ $sbatch_status -ne 0 ]; then
        echo "ERROR: Failed to submit job: $sbatch_output" | tee -a "$MAIN_LOG"
        return 1
    fi
    
    # Extract job ID with robust pattern matching
    local job_id=$(echo "$sbatch_output" | grep -oP 'Submitted batch job \K[0-9]+' || echo "")
    
    if [ -z "$job_id" ]; then
        echo "WARNING: Could not extract job ID from sbatch output: $sbatch_output" | tee -a "$MAIN_LOG"
        return 1
    fi
    
    echo "Submitted ${job_name} job: $job_id" | tee -a "$MAIN_LOG"
    echo "$job_id"
}

# Submit all jobs using the function
echo "Submitting bulk RNA generation job..." | tee -a "$MAIN_LOG"
BULK_JOB_ID=$(submit_job "bulk_gen" "Rscript $SCRIPT_DIR/bulk_generation.R $RLIBRARY $BULK_FILE $SUBDIR_PATH $PREFIX" 8 "16G" "1:00:00")

echo "Submitting single cell labels generation job..." | tee -a "$MAIN_LOG"
LABELS_JOB_ID=$(submit_job "singleCellLabels_gen" "Rscript $SCRIPT_DIR/singleCellLabels_generation.R $RLIBRARY $SEURAT_FILE $SUBDIR_PATH $PREFIX" 8 "16G" "1:00:00")

echo "Submitting single cell expression generation job..." | tee -a "$MAIN_LOG"
SCRNA_JOB_ID=$(submit_job "singleCellExpr_gen" "Rscript $SCRIPT_DIR/singleCellExpr_generation.R $RLIBRARY $SEURAT_FILE $SUBDIR_PATH $MAPPING_FILE $PREFIX" 16 "48G" "2:00:00")

echo "Submitting single cell subjects generation job..." | tee -a "$MAIN_LOG"
SUBJECTS_JOB_ID=$(submit_job "singleCellSubjects_gen" "Rscript $SCRIPT_DIR/singleCellSubjects_generation.R $RLIBRARY $SEURAT_FILE $SUBDIR_PATH $PREFIX" 8 "16G" "1:00:00")

echo "Submitting ground truth generation job..." | tee -a "$MAIN_LOG"b
GT_JOB_ID=$(submit_job "GT_gen" "Rscript $SCRIPT_DIR/GT_generation.R $RLIBRARY $SUBDIR_PATH $SUBDIR_PATH $PREFIX" 8 "16G" "1:00:00")

echo "Submitting per-sample ground truth generation job..." | tee -a "$MAIN_LOG"
SAMPLE_GT_JOB_ID=$(submit_job "GT_per_sample_gen" "Rscript $SCRIPT_DIR/GT_per_sample_generation.R $RLIBRARY $SEURAT_FILE $SUBDIR_PATH $PREFIX" 8 "32G" "1:00:00")

# Create consistent logs directory structure
GLOBAL_LOG_DIR="/work/gr-fe/lorthiois/DeconBenchmark/logs/${PREFIX}"
mkdir -p $GLOBAL_LOG_DIR

# Change where the job mapping file is saved
JOB_MAPPING_FILE="$GLOBAL_LOG_DIR/data_job_mapping.txt"
echo "# Job mapping file for ${PREFIX} data generation - Created $(date)" > "$JOB_MAPPING_FILE"
[ -n "$BULK_JOB_ID" ] && echo "bulk_gen:${BULK_JOB_ID}" >> "$JOB_MAPPING_FILE"
[ -n "$LABELS_JOB_ID" ] && echo "singleCellLabels_gen:${LABELS_JOB_ID}" >> "$JOB_MAPPING_FILE"
[ -n "$SCRNA_JOB_ID" ] && echo "singleCellExpr_gen:${SCRNA_JOB_ID}" >> "$JOB_MAPPING_FILE"
[ -n "$SUBJECTS_JOB_ID" ] && echo "singleCellSubjects_gen:${SUBJECTS_JOB_ID}" >> "$JOB_MAPPING_FILE"
[ -n "$GT_JOB_ID" ] && echo "GT_gen:${GT_JOB_ID}" >> "$JOB_MAPPING_FILE"
[ -n "$SAMPLE_GT_JOB_ID" ] && echo "GT_per_sample_gen:${SAMPLE_GT_JOB_ID}" >> "$JOB_MAPPING_FILE"

echo "All data generation jobs submitted successfully." | tee -a "$MAIN_LOG"
echo "All generated files will be stored in: $SUBDIR_PATH" | tee -a "$MAIN_LOG"
echo "All log files will be stored in: $SCRATCH_LOG_DIR" | tee -a "$MAIN_LOG"
echo "All temporary scripts are stored in: $SCRATCH_SCRIPT_DIR" | tee -a "$MAIN_LOG"
echo "Job mapping saved to: $JOB_MAPPING_FILE" | tee -a "$MAIN_LOG"