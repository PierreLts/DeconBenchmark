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

# Log file setup
LOG_DIR="$SUBDIR_PATH/logs"
mkdir -p $LOG_DIR
MAIN_LOG="$LOG_DIR/master_data_generation_$(date +%Y%m%d_%H%M%S).log"

echo "==== Starting master data generation $(date) ====" | tee -a "$MAIN_LOG"
echo "Seurat file: $SEURAT_FILE" | tee -a "$MAIN_LOG"
echo "Bulk file: $BULK_FILE" | tee -a "$MAIN_LOG" 
echo "Mapping file: $MAPPING_FILE" | tee -a "$MAIN_LOG"
echo "Output directory: $SUBDIR_PATH" | tee -a "$MAIN_LOG"
echo "File prefix: $PREFIX" | tee -a "$MAIN_LOG"
echo "R library path: $RLIBRARY" | tee -a "$MAIN_LOG"
echo "" | tee -a "$MAIN_LOG"

# Base directory where scripts are located
SCRIPT_DIR="/work/gr-fe/lorthiois/DeconBenchmark/scripts/data"

# Submit all jobs independently
echo "Submitting bulk RNA generation job..." | tee -a "$MAIN_LOG"
BULK_JOB_ID=$(sbatch --job-name="bulk_gen_${PREFIX}" \
    --output="$LOG_DIR/%j.out" \
    --error="$LOG_DIR/%j.err" \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task=8 \
    --mem=16G \
    --time=1:00:00 \
    --wrap="module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/; module load r; Rscript $SCRIPT_DIR/bulkRNA_generation.R $RLIBRARY $BULK_FILE $SUBDIR_PATH $PREFIX" \
    | grep -o '[0-9]*')
echo "Submitted bulk job: $BULK_JOB_ID" | tee -a "$MAIN_LOG"

echo "Submitting labels RNA generation job..." | tee -a "$MAIN_LOG"
LABELS_JOB_ID=$(sbatch --job-name="labels_gen_${PREFIX}" \
    --output="$LOG_DIR/%j.out" \
    --error="$LOG_DIR/%j.err" \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task=8 \
    --mem=16G \
    --time=1:00:00 \
    --wrap="module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/; module load r; Rscript $SCRIPT_DIR/labelsRNA_generation.R $RLIBRARY $SEURAT_FILE $SUBDIR_PATH $PREFIX" \
    | grep -o '[0-9]*')
echo "Submitted labels job: $LABELS_JOB_ID" | tee -a "$MAIN_LOG"

echo "Submitting scRNA generation job..." | tee -a "$MAIN_LOG"
SCRNA_JOB_ID=$(sbatch --job-name="scRNA_gen_${PREFIX}" \
    --output="$LOG_DIR/%j.out" \
    --error="$LOG_DIR/%j.err" \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task=16 \
    --mem=48G \
    --time=2:00:00 \
    --wrap="module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/; module load r; Rscript $SCRIPT_DIR/scRNA_generation.R $RLIBRARY $SEURAT_FILE $SUBDIR_PATH $MAPPING_FILE $PREFIX" \
    | grep -o '[0-9]*')
echo "Submitted scRNA job: $SCRNA_JOB_ID" | tee -a "$MAIN_LOG"

echo "Submitting subject information generation job..." | tee -a "$MAIN_LOG"
SUBJECTS_JOB_ID=$(sbatch --job-name="subjects_gen_${PREFIX}" \
    --output="$LOG_DIR/%j.out" \
    --error="$LOG_DIR/%j.err" \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task=8 \
    --mem=16G \
    --time=1:00:00 \
    --wrap="module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/; module load r; Rscript $SCRIPT_DIR/seurat_subjects_generation.R $RLIBRARY $SEURAT_FILE $SUBDIR_PATH $PREFIX" \
    | grep -o '[0-9]*')
echo "Submitted subjects job: $SUBJECTS_JOB_ID" | tee -a "$MAIN_LOG"

echo "Submitting marker genes generation job..." | tee -a "$MAIN_LOG"
MARKER_JOB_ID=$(sbatch --job-name="markers_gen_${PREFIX}" \
    --output="$LOG_DIR/%j.out" \
    --error="$LOG_DIR/%j.err" \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task=16 \
    --mem=32G \
    --time=2:00:00 \
    --wrap="module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/; module load r; Rscript $SCRIPT_DIR/seurat_markers_generation.R $RLIBRARY $SEURAT_FILE $SUBDIR_PATH $MAPPING_FILE $PREFIX" \
    | grep -o '[0-9]*')
echo "Submitted markers job: $MARKER_JOB_ID" | tee -a "$MAIN_LOG"

# Only submit dependent job if we have a valid MARKER_JOB_ID
if [[ -n "$MARKER_JOB_ID" && "$MARKER_JOB_ID" =~ ^[0-9]+$ ]]; then
    echo "Submitting significant genes generation job (dependent on markers)..." | tee -a "$MAIN_LOG"
    SIG_GENES_JOB_ID=$(sbatch --job-name="sig_genes_gen_${PREFIX}" \
        --output="$LOG_DIR/%j.out" \
        --error="$LOG_DIR/%j.err" \
        --nodes=1 \
        --ntasks=1 \
        --cpus-per-task=8 \
        --mem=16G \
        --time=1:00:00 \
        --dependency=afterok:$MARKER_JOB_ID \
        --wrap="module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/; module load r; Rscript $SCRIPT_DIR/seurat_significant_genes_generation.R $RLIBRARY $SUBDIR_PATH/${PREFIX}_markers.rda $SUBDIR_PATH $PREFIX" \
        | grep -o '[0-9]*')
    echo "Submitted significant genes job: $SIG_GENES_JOB_ID" | tee -a "$MAIN_LOG"
else
    echo "WARNING: Invalid marker job ID. Skipping significant genes generation." | tee -a "$MAIN_LOG"
fi

echo "Submitting cell type expression generation job..." | tee -a "$MAIN_LOG"
CELLTYPE_JOB_ID=$(sbatch --job-name="celltype_expr_gen_${PREFIX}" \
    --output="$LOG_DIR/%j.out" \
    --error="$LOG_DIR/%j.err" \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task=16 \
    --mem=32G \
    --time=2:00:00 \
    --wrap="module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/; module load r; Rscript $SCRIPT_DIR/seurat_celltype_expression_generation.R $RLIBRARY $SEURAT_FILE $MAPPING_FILE $SUBDIR_PATH $PREFIX" \
    | grep -o '[0-9]*')
echo "Submitted cell type expression job: $CELLTYPE_JOB_ID" | tee -a "$MAIN_LOG"

echo "Submitting ground truth generation job..." | tee -a "$MAIN_LOG"
GT_JOB_ID=$(sbatch --job-name="ground_truth_${PREFIX}" \
    --output="$LOG_DIR/%j.out" \
    --error="$LOG_DIR/%j.err" \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task=8 \
    --mem=16G \
    --time=1:00:00 \
    --wrap="module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/; module load r; Rscript $SCRIPT_DIR/ground_truth.R $RLIBRARY $SUBDIR_PATH $SUBDIR_PATH $PREFIX" \
    | grep -o '[0-9]*')
echo "Submitted ground truth job: $GT_JOB_ID" | tee -a "$MAIN_LOG"

echo "Submitting per-sample ground truth generation job..." | tee -a "$MAIN_LOG"
SAMPLE_GT_JOB_ID=$(sbatch --job-name="sample_gt_gen_${PREFIX}" \
    --output="$LOG_DIR/%j.out" \
    --error="$LOG_DIR/%j.err" \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task=8 \
    --mem=32G \
    --time=1:00:00 \
    --wrap="module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/; module load r; Rscript $SCRIPT_DIR/per_sample_ground_truth.R $RLIBRARY $SEURAT_FILE $SUBDIR_PATH $PREFIX" \
    | grep -o '[0-9]*')
echo "Submitted per-sample ground truth job: $SAMPLE_GT_JOB_ID" | tee -a "$MAIN_LOG"

echo "All data generation jobs submitted successfully." | tee -a "$MAIN_LOG"
echo "All generated files will be stored in: $SUBDIR_PATH" | tee -a "$MAIN_LOG"