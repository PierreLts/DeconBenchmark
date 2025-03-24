#!/bin/bash
# parallel_deconv.sh - Run multiple deconvolution methods in parallel

# Default parameters
DEFAULT_DATASET_PREFIX="TB_D33348"
DEFAULT_SAMPLE_FILTER="AB"
DEFAULT_BULK_TYPE="bulk"  # bulk, bulk_random, pseudobulk
DEFAULT_METHODS="ARIC,AutoGeneS,BayesPrism,BisqueRef,DeconPeaker,DeconRNASeq,DESeq2,EMeth,EPIC,FARDEEP,LinDeconSeq,MIXTURE,MuSic,MySort,PREDE,quanTIseq,RNA-Sieve,scaden,SCDC,TOAST"

# Working ref, but long: "deconvSeq,DWLS,CPM,BayICE"
# Working ref free: "BayCount,CDSeq,deconf,DeconICA,DSA,Linseed,MCPcounter"
# Entire list: "AdRoit,ARIC,AutoGeneS,BayCount,BayesPrism,BayICE,BisqueMarker,BisqueRef,BseqSC,CDSeq,CellDistinguisher,CIBERSORT,CIBERSORTx,CPM,DAISM,debCAM,Deblender,DeCompress,deconf,DeconICA,DeconPeaker,DeconRNASeq,deconvSeq,DecOT,DeMixT,DESeq2,digitalDLSorter,DSA,dtangle,DWLS,EMeth,EPIC,FARDEEP,ImmuCellAI,LinDeconSeq,Linseed,MCPcounter,MIXTURE,MOMF,MuSic,MySort,NITUMID,PREDE,quanTIseq,RNA-Sieve,scaden,SCDC,spatialDWLS,TOAST"

# Parse command line arguments
DATASET_PREFIX="${1:-$DEFAULT_DATASET_PREFIX}"
SAMPLE_FILTER="${2:-$DEFAULT_SAMPLE_FILTER}"
METHODS="${3:-$DEFAULT_METHODS}"
BULK_TYPE="${4:-$DEFAULT_BULK_TYPE}"  # New parameter - bulk file type

# Create the new output directory name with dataset prefix and bulk type
OUTPUT_DIR_NAME="${DATASET_PREFIX}-${BULK_TYPE}"

# Set paths
SCRIPT_DIR="/scratch/lorthiois/scripts"
mkdir -p "$SCRIPT_DIR"
TEMPLATE_DIR="/work/gr-fe/lorthiois/DeconBenchmark/scripts/deconvolution"
OUTPUT_BASE_DIR="/work/gr-fe/lorthiois/DeconBenchmark/deconv_results"
OUTPUT_DIR="$OUTPUT_BASE_DIR/$OUTPUT_DIR_NAME"
BENCHMARK_DIR="/work/gr-fe/lorthiois/DeconBenchmark/benchmark_results/$OUTPUT_DIR_NAME"
GLOBAL_LOG_DIR="/work/gr-fe/lorthiois/DeconBenchmark/logs/${OUTPUT_DIR_NAME}_${SAMPLE_FILTER}"
mkdir -p ${GLOBAL_LOG_DIR}

# Create output directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$BENCHMARK_DIR"

# Log file setup
LOG_DIR="$OUTPUT_DIR/logs"
mkdir -p "$LOG_DIR"
MAIN_LOG="$LOG_DIR/parallel_deconv_${SAMPLE_FILTER}_$(date +%Y%m%d_%H%M%S).log"

echo "==== Starting parallel deconvolution pipeline $(date) ====" | tee -a "$MAIN_LOG"
echo "Dataset prefix: $DATASET_PREFIX" | tee -a "$MAIN_LOG"
echo "Sample filter: $SAMPLE_FILTER" | tee -a "$MAIN_LOG"
echo "Bulk file type: $BULK_TYPE" | tee -a "$MAIN_LOG"
echo "Output directory name: $OUTPUT_DIR_NAME" | tee -a "$MAIN_LOG"
echo "Methods to run: ${METHODS}" | tee -a "$MAIN_LOG"
echo "Results will be saved to: $OUTPUT_DIR" | tee -a "$MAIN_LOG"
echo "" | tee -a "$MAIN_LOG"
echo "Storing model job mapping in: $GLOBAL_LOG_DIR/model_job_mapping.txt" | tee -a "$MAIN_LOG"

# Convert methods string to array
IFS=',' read -r -a MODELS <<< "$METHODS"
echo "Total methods: ${#MODELS[@]}" | tee -a "$MAIN_LOG"

# Array to store job IDs for each model's deconvolution job
declare -a DECONV_JOB_IDS

# Submit deconvolution jobs for all methods in parallel
for MODEL in "${MODELS[@]}"; do
    # Skip CDSeq in the main loop since it will be handled separately
    if [[ "$MODEL" == "CDSeq" ]]; then
        echo "Skipping CDSeq in main loop - will be run with special script" | tee -a "$MAIN_LOG"
        continue
    fi
    
    echo "Submitting deconvolution job for model: $MODEL" | tee -a "$MAIN_LOG"
    
    # Create a temporary job script for this model
    TEMP_SCRIPT="$SCRIPT_DIR/temp_${MODEL}_${OUTPUT_DIR_NAME}_${SAMPLE_FILTER}_deconv.sh"
    
    # Copy the template and modify it to directly use our output directory
    sed -e "s|#SBATCH --job-name=.*|#SBATCH --job-name=${MODEL}_${OUTPUT_DIR_NAME}_${SAMPLE_FILTER}|" \
        -e "s|#SBATCH --time=.*|#SBATCH --time=48:00:00|" \
        -e "s|dataset_prefix=.*|dataset_prefix=\"$DATASET_PREFIX\"|" \
        -e "s|sample_filter=.*|sample_filter=\"$SAMPLE_FILTER\"|" \
        -e "s|deconv_method=.*|deconv_method=\"$MODEL\"|" \
        -e "s|bulk_type=.*|bulk_type=\"$BULK_TYPE\"|" \
        -e "s|output_base_dir=.*|output_base_dir=\"$OUTPUT_BASE_DIR\"|" \
        -e "s|output_dir=.*|output_dir=\"$OUTPUT_DIR\"|" \
        "$TEMPLATE_DIR/deconv_run.sh" > "$TEMP_SCRIPT"
    
    # Add an override at the end of the script to ensure results go to the right place
    cat >> "$TEMP_SCRIPT" << EOF

# Override Rscript call to force output directory
DEST_DIR="$OUTPUT_DIR"
Rscript \${SCRIPT} \${RLIBRARY} \${dataset_prefix} \${sample_filter} \${DEST_DIR} \${deconv_method} \${sparse_conversion} \${bulk_type}
exit 0  # Exit after our custom Rscript call
EOF
    
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
    echo "${MODEL}:${JOB_ID}" >> "$GLOBAL_LOG_DIR/model_job_mapping.txt"
done

# Special case for CDSeq which requires a separate execution
if [[ "$METHODS" == *"CDSeq"* ]] || [[ "$METHODS" == *"All"* ]]; then
    echo "Submitting special CDSeq job" | tee -a "$MAIN_LOG"
    CDSEQ_SCRIPT="$SCRIPT_DIR/temp_cdseq_${OUTPUT_DIR_NAME}_${SAMPLE_FILTER}_deconv.sh"
    
    # Create CDSeq script by modifying the template and adding our output directory override
    sed -e "s|#SBATCH --job-name=.*|#SBATCH --job-name=CDSeq_${OUTPUT_DIR_NAME}_${SAMPLE_FILTER}|" \
        -e "s|dataset_prefix=.*|dataset_prefix=\"$DATASET_PREFIX\"|" \
        -e "s|sample_filter=.*|sample_filter=\"$SAMPLE_FILTER\"|" \
        -e "s|bulk_type=.*|bulk_type=\"$BULK_TYPE\"|" \
        -e "s|output_base_dir=.*|output_base_dir=\"$OUTPUT_BASE_DIR\"|" \
        -e "s|output_dir=.*|output_dir=\"$OUTPUT_DIR\"|" \
        "$TEMPLATE_DIR/run_cdseq.sh" > "$CDSEQ_SCRIPT"
    
    # Add override at the end
    cat >> "$CDSEQ_SCRIPT" << EOF

# Override Rscript call to force output directory
DEST_DIR="$OUTPUT_DIR"
Rscript \${SCRIPT} \${RLIBRARY} \${input_data} \${DEST_DIR} \${sample_filter} \${bulk_type}
exit 0  # Exit after our custom Rscript call
EOF
    
    chmod +x "$CDSEQ_SCRIPT"
    
    # Submit CDSeq job
    CDSEQ_JOB_ID=$(sbatch "$CDSEQ_SCRIPT" | grep -o '[0-9]*')
    if [ -n "$CDSEQ_JOB_ID" ]; then
        echo "Submitted CDSeq job $CDSEQ_JOB_ID" | tee -a "$MAIN_LOG"
        DECONV_JOB_IDS+=("$CDSEQ_JOB_ID")
        echo "CDSeq:${CDSEQ_JOB_ID}" >> "$GLOBAL_LOG_DIR/model_job_mapping.txt"
    else
        echo "ERROR: Failed to submit CDSeq job" | tee -a "$MAIN_LOG"
    fi
fi

echo "All deconvolution jobs submitted successfully." | tee -a "$MAIN_LOG"
echo "Results will be available in: $OUTPUT_DIR" | tee -a "$MAIN_LOG"