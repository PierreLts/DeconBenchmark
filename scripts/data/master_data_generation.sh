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

# Create prefix-specific subdirectory
SUBDIR_PATH="$OUTPUT_DIR/$PREFIX"
mkdir -p $SUBDIR_PATH

# Set up log directories
SCRATCH_LOG_DIR="/scratch/lorthiois/logs"
mkdir -p $SCRATCH_LOG_DIR
LOCAL_LOG_DIR="$SUBDIR_PATH/logs"
mkdir -p $LOCAL_LOG_DIR
MAIN_LOG="$LOCAL_LOG_DIR/master_data_generation_$(date +%Y%m%d_%H%M%S).log"

echo "==== Starting master data generation $(date) ====" | tee -a "$MAIN_LOG"
echo "Seurat file: $SEURAT_FILE" | tee -a "$MAIN_LOG"
echo "Bulk file: $BULK_FILE" | tee -a "$MAIN_LOG" 
echo "Mapping file: $MAPPING_FILE" | tee -a "$MAIN_LOG"
echo "Output directory: $SUBDIR_PATH" | tee -a "$MAIN_LOG"
echo "File prefix: $PREFIX" | tee -a "$MAIN_LOG"
echo "R library path: $RLIBRARY" | tee -a "$MAIN_LOG"
echo "Log directory: $SCRATCH_LOG_DIR" | tee -a "$MAIN_LOG"
echo "" | tee -a "$MAIN_LOG"

# Base directory where scripts are located
SCRIPT_DIR="/work/gr-fe/lorthiois/DeconBenchmark/scripts/data"

# Function to submit a job
submit_job() {
    local script=$1
    local job_name=$2
    local args="${@:3}"  # All arguments after the first two
    
    echo "Submitting job: $job_name" | tee -a "$MAIN_LOG"
    echo "Command: Rscript $script $args" | tee -a "$MAIN_LOG"
    
    # Submit job using sbatch with logs redirected to scratch directory
    JOB_ID=$(sbatch --job-name=$job_name \
                    --output=$SCRATCH_LOG_DIR/%j.o \
                    --error=$SCRATCH_LOG_DIR/%j.e \
                    --nodes=1 \
                    --ntasks=1 \
                    --cpus-per-task=8 \
                    --mem=32G \
                    --time=2:00:00 \
                    --wrap="module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/; module load r; Rscript $script $args" \
            | grep -o '[0-9]*')
    
    if [ -z "$JOB_ID" ]; then
        echo "ERROR: Failed to submit job $job_name" | tee -a "$MAIN_LOG"
    else
        echo "Submitted job $JOB_ID for $job_name" | tee -a "$MAIN_LOG"
    fi
    
    echo $JOB_ID
}

# Submit all data generation jobs in parallel
JOBS=()

# 1. bulkRNA_generation
JOB_ID=$(submit_job "$SCRIPT_DIR/bulkRNA_generation.R" \
                   "bulk_gen" \
                   "$RLIBRARY" \
                   "$BULK_FILE" \
                   "$SUBDIR_PATH" \
                   "$PREFIX")
JOBS+=($JOB_ID)

# 2. labelsRNA_generation
JOB_ID=$(submit_job "$SCRIPT_DIR/labelsRNA_generation.R" \
                   "labels_gen" \
                   "$RLIBRARY" \
                   "$SEURAT_FILE" \
                   "$SUBDIR_PATH" \
                   "$PREFIX")
JOBS+=($JOB_ID)

# 3. scRNA_generation
JOB_ID=$(submit_job "$SCRIPT_DIR/scRNA_generation.R" \
                   "scRNA_gen" \
                   "$RLIBRARY" \
                   "$SEURAT_FILE" \
                   "$SUBDIR_PATH" \
                   "$MAPPING_FILE" \
                   "$PREFIX")
JOBS+=($JOB_ID)

# 4. Subject information generation
JOB_ID=$(submit_job "$SCRIPT_DIR/seurat_subjects_generation.R" \
                   "subjects_gen" \
                   "$RLIBRARY" \
                   "$SEURAT_FILE" \
                   "$SUBDIR_PATH" \
                   "$PREFIX")
JOBS+=($JOB_ID)

# 5. Marker genes generation
JOB_ID=$(submit_job "$SCRIPT_DIR/seurat_markers_generation.R" \
                   "markers_gen" \
                   "$RLIBRARY" \
                   "$SEURAT_FILE" \
                   "$SUBDIR_PATH" \
                   "$MAPPING_FILE" \
                   "$PREFIX")
MARKER_JOB_ID=$JOB_ID
JOBS+=($JOB_ID)

# 6. Wait for marker job to complete before running significant genes
if [ ! -z "$MARKER_JOB_ID" ]; then
    # Create a temporary script that waits for the marker job
    TEMP_SCRIPT="$LOCAL_LOG_DIR/run_significant_genes_after_markers.sh"
    cat > $TEMP_SCRIPT << EOF
#!/bin/bash
#SBATCH --job-name=sig_genes_gen
#SBATCH --output=$SCRATCH_LOG_DIR/%j.o
#SBATCH --error=$SCRATCH_LOG_DIR/%j.e
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --dependency=afterok:$MARKER_JOB_ID

module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/
module load r
Rscript $SCRIPT_DIR/seurat_significant_genes_generation.R $RLIBRARY "$SUBDIR_PATH/${PREFIX}_markers.rda" $SUBDIR_PATH $PREFIX
EOF

    chmod +x $TEMP_SCRIPT
    SIG_GENES_JOB_ID=$(sbatch $TEMP_SCRIPT | grep -o '[0-9]*')
    echo "Submitted job $SIG_GENES_JOB_ID for significant genes generation (dependent on markers job)" | tee -a "$MAIN_LOG"
    JOBS+=($SIG_GENES_JOB_ID)
fi

# 7. Cell type expression generation
JOB_ID=$(submit_job "$SCRIPT_DIR/seurat_celltype_expression_generation.R" \
                   "celltype_expr_gen" \
                   "$RLIBRARY" \
                   "$SEURAT_FILE" \
                   "$MAPPING_FILE" \
                   "$SUBDIR_PATH" \
                   "$PREFIX")
JOBS+=($JOB_ID)

# 8. Ground truth generation - depends on single cell labels being ready
TEMP_SCRIPT="$LOCAL_LOG_DIR/run_ground_truth.sh"
cat > $TEMP_SCRIPT << EOF
#!/bin/bash
#SBATCH --job-name=ground_truth
#SBATCH --output=$SCRATCH_LOG_DIR/%j.o
#SBATCH --error=$SCRATCH_LOG_DIR/%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1:00:00

module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/
module load r
Rscript $SCRIPT_DIR/ground_truth.R $RLIBRARY $SUBDIR_PATH $SUBDIR_PATH $PREFIX
EOF

chmod +x $TEMP_SCRIPT
GT_JOB_ID=$(sbatch $TEMP_SCRIPT | grep -o '[0-9]*')
echo "Submitted job $GT_JOB_ID for ground truth generation" | tee -a "$MAIN_LOG"
JOBS+=($GT_JOB_ID)

# 9. Per-sample ground truth generation
JOB_ID=$(submit_job "$SCRIPT_DIR/per_sample_ground_truth.R" \
                   "sample_gt_gen" \
                   "$RLIBRARY" \
                   "$SEURAT_FILE" \
                   "$SUBDIR_PATH" \
                   "$PREFIX")
JOBS+=($JOB_ID)

# 10. Final RDA generation - depends on bulk, single-cell expression, and labels being ready
# Create comma-separated list of job IDs for dependency
JOB_DEPS=$(IFS=,; echo "${JOBS[*]}")

FINAL_SCRIPT="$LOCAL_LOG_DIR/run_final_rda.sh"
cat > $FINAL_SCRIPT << EOF
#!/bin/bash
#SBATCH --job-name=finalRDA
#SBATCH --output=$SCRATCH_LOG_DIR/%j.o
#SBATCH --error=$SCRATCH_LOG_DIR/%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --dependency=afterok:$JOB_DEPS

module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/
module load r
Rscript $SCRIPT_DIR/finalRDA_generation.R $RLIBRARY $SUBDIR_PATH $SUBDIR_PATH $PREFIX
EOF

chmod +x $FINAL_SCRIPT
FINAL_JOB_ID=$(sbatch $FINAL_SCRIPT | grep -o '[0-9]*')
echo "Submitted job $FINAL_JOB_ID for final RDA generation (dependent on all previous jobs)" | tee -a "$MAIN_LOG"

echo "All data generation jobs submitted. Check individual logs in $SCRATCH_LOG_DIR for progress." | tee -a "$MAIN_LOG"
echo "Final RDA will be generated after all other jobs complete." | tee -a "$MAIN_LOG"
echo "All generated files will be stored in: $SUBDIR_PATH" | tee -a "$MAIN_LOG"