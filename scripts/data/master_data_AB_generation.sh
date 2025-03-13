#!/bin/bash
# master_data_AB_generation.sh - Run data generation tasks with A/B/AB filtering

# ===== EXPLICIT SETTINGS - MODIFY THESE DIRECTLY =====
# Output prefix for generated files
PREFIX="TB1"
# Sample filter: "A" for only A samples, "B" for only B samples, "AB" for both
SAMPLE_FILTER="AB"
# =====

# Default parameters
DEFAULT_SEURAT_FILE="/work/gr-fe/lorthiois/DeconBenchmark/data/GFB-33245_HFKJMDSXC_2_scRNAseqWTATBseverityrun1_Seurat.rds"
DEFAULT_BULK_FILE="/work/gr-fe/lorthiois/DeconBenchmark/data/cleaned_feature_counts_matrix.csv"
DEFAULT_MAPPING_FILE="/work/gr-fe/lorthiois/DeconBenchmark/data/mart_export.txt"
DEFAULT_OUTPUT_DIR="/work/gr-fe/lorthiois/DeconBenchmark/generated_data"
DEFAULT_RLIBRARY="/work/gr-fe/R_4.3.1"

# Parse command line arguments
SEURAT_FILE="${1:-$DEFAULT_SEURAT_FILE}"
BULK_FILE="${2:-$DEFAULT_BULK_FILE}"
MAPPING_FILE="${3:-$DEFAULT_MAPPING_FILE}"
OUTPUT_DIR="${4:-$DEFAULT_OUTPUT_DIR}"
# Override PREFIX if provided as 5th argument
if [ -n "$5" ]; then
    PREFIX="$5"
fi
# Override SAMPLE_FILTER if provided as 6th argument
if [ -n "$6" ]; then
    SAMPLE_FILTER="$6"
fi
RLIBRARY="${7:-$DEFAULT_RLIBRARY}"

# Output prefix for filtered data
FILTERED_PREFIX="${PREFIX}-${SAMPLE_FILTER}"

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
MAIN_LOG="$OUTPUT_LOG_DIR/master_data_AB_generation_$(date +%Y%m%d_%H%M%S).log"

echo "==== Starting master data AB generation $(date) ====" | tee -a "$MAIN_LOG"
echo "Seurat file: $SEURAT_FILE" | tee -a "$MAIN_LOG"
echo "Bulk file: $BULK_FILE" | tee -a "$MAIN_LOG" 
echo "Mapping file: $MAPPING_FILE" | tee -a "$MAIN_LOG"
echo "Output directory: $SUBDIR_PATH" | tee -a "$MAIN_LOG"
echo "File prefix: $PREFIX" | tee -a "$MAIN_LOG"
echo "Sample filter: $SAMPLE_FILTER" | tee -a "$MAIN_LOG"
echo "Filtered prefix: $FILTERED_PREFIX" | tee -a "$MAIN_LOG"
echo "R library path: $RLIBRARY" | tee -a "$MAIN_LOG"
echo "Scratch log directory: $SCRATCH_LOG_DIR" | tee -a "$MAIN_LOG"
echo "Scratch script directory: $SCRATCH_SCRIPT_DIR" | tee -a "$MAIN_LOG"
echo "" | tee -a "$MAIN_LOG"

# Validate filter value
if [[ ! "$SAMPLE_FILTER" =~ ^(A|B|AB)$ ]]; then
    echo "ERROR: SAMPLE_FILTER must be one of 'A', 'B', or 'AB'. Got '$SAMPLE_FILTER'" | tee -a "$MAIN_LOG"
    exit 1
fi

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
    local temp_script="$SCRATCH_SCRIPT_DIR/${job_name}_${PREFIX}_${SAMPLE_FILTER}.sh"
    echo "#!/bin/bash" > "$temp_script"
    echo "#SBATCH --job-name=${job_name}_${PREFIX}_${SAMPLE_FILTER}" >> "$temp_script"
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
echo "Submitting bulk RNA generation job (unfiltered)..." | tee -a "$MAIN_LOG"
BULK_JOB_ID=$(submit_job "bulk_gen" "Rscript $SCRIPT_DIR/bulk_generation.R $RLIBRARY $BULK_FILE $SUBDIR_PATH $PREFIX" 8 "16G" "1:00:00")

echo "Submitting single cell labels generation job with $SAMPLE_FILTER filter..." | tee -a "$MAIN_LOG"
LABELS_JOB_ID=$(submit_job "singleCellLabels_gen" "
# Explicit settings for AB generation
PREFIX=\"$PREFIX\"
SAMPLE_FILTER=\"$SAMPLE_FILTER\"

# Execute the filtered version
Rscript $SCRIPT_DIR/singleCellLabels_AB_generation.R $RLIBRARY $SEURAT_FILE $SUBDIR_PATH $PREFIX $SAMPLE_FILTER
" 8 "16G" "1:00:00")

echo "Submitting single cell expression generation job with $SAMPLE_FILTER filter..." | tee -a "$MAIN_LOG"
SCRNA_JOB_ID=$(submit_job "singleCellExpr_gen" "
# Explicit settings for AB generation
PREFIX=\"$PREFIX\"
SAMPLE_FILTER=\"$SAMPLE_FILTER\"

# Execute the filtered version
Rscript $SCRIPT_DIR/singleCellExpr_AB_generation.R $RLIBRARY $SEURAT_FILE $SUBDIR_PATH $MAPPING_FILE $PREFIX $SAMPLE_FILTER
" 16 "48G" "2:00:00")

echo "Submitting single cell subjects generation job with $SAMPLE_FILTER filter..." | tee -a "$MAIN_LOG"
SUBJECTS_JOB_ID=$(submit_job "singleCellSubjects_gen" "
# Explicit settings for AB generation
PREFIX=\"$PREFIX\"
SAMPLE_FILTER=\"$SAMPLE_FILTER\"

# Execute the filtered version
Rscript $SCRIPT_DIR/singleCellSubjects_AB_generation.R $RLIBRARY $SEURAT_FILE $SUBDIR_PATH $PREFIX $SAMPLE_FILTER
" 8 "16G" "1:00:00")


echo "Submitting pseudobulk transfer job..." | tee -a "$MAIN_LOG"
PSEUDOBULK_JOB_ID=$(submit_job "transfer_pseudobulk" "
# Process pseudobulk data with gene name mapping
Rscript $SCRIPT_DIR/transfer_pseudobulk.R $RLIBRARY /work/gr-fe/lorthiois/DeconBenchmark/data/pseudobulk_counts_120k.csv $MAPPING_FILE $OUTPUT_BASE_DIR $PREFIX $SAMPLE_FILTER
" 4 "8G" "0:30:00")


# Create consistent logs directory structure
GLOBAL_LOG_DIR="/work/gr-fe/lorthiois/DeconBenchmark/logs/${PREFIX}_${SAMPLE_FILTER}"
mkdir -p $GLOBAL_LOG_DIR

# Change where the job mapping file is saved
JOB_MAPPING_FILE="$GLOBAL_LOG_DIR/data_job_mapping.txt"
echo "# Job mapping file for ${PREFIX}_${SAMPLE_FILTER} data generation - Created $(date)" > "$JOB_MAPPING_FILE"
[ -n "$BULK_JOB_ID" ] && echo "bulk_gen:${BULK_JOB_ID}" >> "$JOB_MAPPING_FILE"
[ -n "$LABELS_JOB_ID" ] && echo "singleCellLabels_gen:${LABELS_JOB_ID}" >> "$JOB_MAPPING_FILE"
[ -n "$SCRNA_JOB_ID" ] && echo "singleCellExpr_gen:${SCRNA_JOB_ID}" >> "$JOB_MAPPING_FILE"
[ -n "$SUBJECTS_JOB_ID" ] && echo "singleCellSubjects_gen:${SUBJECTS_JOB_ID}" >> "$JOB_MAPPING_FILE"
[ -n "$PSEUDOBULK_JOB_ID" ] && echo "transfer_pseudobulk:${PSEUDOBULK_JOB_ID}" >> "$JOB_MAPPING_FILE"


# Add final job to create a consolidated RDA file when all others have completed
JOBS_DEPENDENCY="${BULK_JOB_ID},${LABELS_JOB_ID},${SCRNA_JOB_ID},${SUBJECTS_JOB_ID},${PSEUDOBULK_JOB_ID}"

echo "Submitting finalRDA generation job..." | tee -a "$MAIN_LOG"
FINAL_RDA_JOB_ID=$(submit_job "finalRDA_gen" "
# Create a single consolidated RDA file for the filtered dataset
Rscript $SCRIPT_DIR/finalRDA_generation.R $RLIBRARY $SUBDIR_PATH $FILTERED_PREFIX

# Print confirmation
echo \"Final RDA generation completed for $FILTERED_PREFIX\"
" 8 "32G" "0:30:00" "$JOBS_DEPENDENCY")

[ -n "$FINAL_RDA_JOB_ID" ] && echo "finalRDA_gen:${FINAL_RDA_JOB_ID}" >> "$JOB_MAPPING_FILE"

echo "All data generation jobs submitted successfully." | tee -a "$MAIN_LOG"
echo "All generated files will be stored in: $SUBDIR_PATH" | tee -a "$MAIN_LOG"
echo "All log files will be stored in: $SCRATCH_LOG_DIR" | tee -a "$MAIN_LOG"
echo "All temporary scripts are stored in: $SCRATCH_SCRIPT_DIR" | tee -a "$MAIN_LOG"
echo "Job mapping saved to: $JOB_MAPPING_FILE" | tee -a "$MAIN_LOG"