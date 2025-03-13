#!/bin/bash
# parallel_deconv.sh - Run multiple deconvolution methods in parallel

# Default parameters
DEFAULT_DATASET_PREFIX="TB1-random"
DEFAULT_SAMPLE_FILTER="AB"
DEFAULT_METHODS="AdRoit,ARIC,AutoGeneS,BayCount,BayesPrism,BayICE,BisqueMarker,BisqueRef,BseqSC,CDSeq,CellDistinguisher,CIBERSORT,CIBERSORTx,CPM,DAISM,debCAM,Deblender,DeCompress,deconf,DeconICA,DeconPeaker,DeconRNASeq,deconvSeq,DecOT,DeMixT,DESeq2,digitalDLSorter,DSA,dtangle,DWLS,EMeth,EPIC,FARDEEP,ImmuCellAI,LinDeconSeq,Linseed,MCPcounter,MIXTURE,MOMF,MuSic,MySort,NITUMID,PREDE,quanTIseq,RNA-Sieve,scaden,SCDC,spatialDWLS,TOAST"
DEFAULT_BULK_FILE="bulk_random" # Default bulk file name (without prefix_ and .rda) "bulk", "bulk_random", etc.

# Parse command line arguments
DATASET_PREFIX="${1:-$DEFAULT_DATASET_PREFIX}"
SAMPLE_FILTER="${2:-$DEFAULT_SAMPLE_FILTER}"
METHODS="${3:-$DEFAULT_METHODS}"
BULK_FILE="${4:-$DEFAULT_BULK_FILE}" # New parameter for bulk file

# Set paths
SCRIPT_DIR="/scratch/lorthiois/scripts"
mkdir -p "$SCRIPT_DIR"
TEMPLATE_DIR="/work/gr-fe/lorthiois/DeconBenchmark/scripts/deconvolution"
OUTPUT_BASE_DIR="/work/gr-fe/lorthiois/DeconBenchmark/deconv_results"
OUTPUT_DIR="$OUTPUT_BASE_DIR/$DATASET_PREFIX"
GLOBAL_LOG_DIR="/work/gr-fe/lorthiois/DeconBenchmark/logs/${DATASET_PREFIX}_${SAMPLE_FILTER}"
mkdir -p ${GLOBAL_LOG_DIR}

# The full bulk file path
BULK_FILE_PATH="/work/gr-fe/lorthiois/DeconBenchmark/generated_data/${DATASET_PREFIX}/${DATASET_PREFIX}_${BULK_FILE}.rda"

# Create output directories
mkdir -p "$OUTPUT_DIR"

# Log file setup
LOG_DIR="$OUTPUT_DIR/logs"
mkdir -p "$LOG_DIR"
MAIN_LOG="$LOG_DIR/parallel_deconv_${SAMPLE_FILTER}_$(date +%Y%m%d_%H%M%S).log"

echo "==== Starting parallel deconvolution pipeline $(date) ====" | tee -a "$MAIN_LOG"
echo "Dataset prefix: $DATASET_PREFIX" | tee -a "$MAIN_LOG"
echo "Sample filter: $SAMPLE_FILTER" | tee -a "$MAIN_LOG"
echo "Bulk file: $BULK_FILE_PATH" | tee -a "$MAIN_LOG"
echo "Methods to run: ${METHODS}" | tee -a "$MAIN_LOG"
echo "Results will be saved to: $OUTPUT_DIR" | tee -a "$MAIN_LOG"
echo "" | tee -a "$MAIN_LOG"

# Check if bulk file exists
if [ ! -f "$BULK_FILE_PATH" ]; then
    echo "ERROR: Bulk file not found at $BULK_FILE_PATH" | tee -a "$MAIN_LOG"
    exit 1
fi

echo "Storing model job mapping in: $GLOBAL_LOG_DIR/model_job_mapping.txt" | tee -a "$MAIN_LOG"

# Convert methods string to array
IFS=',' read -r -a MODELS <<< "$METHODS"
echo "Total methods: ${#MODELS[@]}" | tee -a "$MAIN_LOG"

# Create/clean the model-job mapping file
> "$GLOBAL_LOG_DIR/model_job_mapping.txt"

# Submit deconvolution jobs for all methods in parallel
for MODEL in "${MODELS[@]}"; do
    # Skip CDSeq in the main loop since it will be handled separately
    if [[ "$MODEL" == "CDSeq" ]]; then
        echo "Skipping CDSeq in main loop - will be run with special script" | tee -a "$MAIN_LOG"
        continue
    fi
    
    echo "Submitting deconvolution job for model: $MODEL" | tee -a "$MAIN_LOG"
    
    # Create a temporary job script for this model
    TEMP_SCRIPT="$SCRIPT_DIR/temp_${MODEL}_${DATASET_PREFIX}_${SAMPLE_FILTER}_deconv.sh"
    
    # Copy the template and modify parameters
    cat "$TEMPLATE_DIR/deconv_run.sh" | \
        sed "s|dataset_prefix=.*|dataset_prefix=\"$DATASET_PREFIX\"|" | \
        sed "s|sample_filter=.*|sample_filter=\"$SAMPLE_FILTER\"|" | \
        sed "s|deconv_method=.*|deconv_method=\"$MODEL\"|" | \
        sed "s|#SBATCH --job-name=.*|#SBATCH --job-name=${MODEL}_${DATASET_PREFIX}_${SAMPLE_FILTER}|" | \
        sed "s|#SBATCH --time=.*|#SBATCH --time=48:00:00|" > "$TEMP_SCRIPT"
    
    # Add custom bulk file parameter to the script
    sed -i "s|SCRIPT=.*|SCRIPT=\${SCRIPT} \${RLIBRARY} \${dataset_prefix} \${sample_filter} \${output_base_dir} \${deconv_method} \${sparse_conversion} \"${BULK_FILE}\"|" "$TEMP_SCRIPT"
        
    # Make executable
    chmod +x "$TEMP_SCRIPT"
    
    # Submit job
    JOB_ID=$(sbatch "$TEMP_SCRIPT" | grep -o '[0-9]*')
    if [ -z "$JOB_ID" ]; then
        echo "ERROR: Failed to submit deconvolution job for $MODEL" | tee -a "$MAIN_LOG"
        continue
    fi
    
    echo "Submitted job $JOB_ID for $MODEL deconvolution" | tee -a "$MAIN_LOG"
    
    # Store the model-job ID mapping for later use
    echo "${MODEL}:${JOB_ID}" >> "$GLOBAL_LOG_DIR/model_job_mapping.txt"
done

# Special case for CDSeq which requires a separate execution
if [[ "$METHODS" == *"CDSeq"* ]] || [[ "$METHODS" == *"All"* ]]; then
    echo "Submitting special CDSeq job" | tee -a "$MAIN_LOG"
    CDSEQ_SCRIPT="$SCRIPT_DIR/temp_cdseq_${DATASET_PREFIX}_${SAMPLE_FILTER}_deconv.sh"
    
    # Create CDSeq script by modifying the template
    cat "$TEMPLATE_DIR/run_cdseq.sh" | \
        sed "s|dataset_prefix=.*|dataset_prefix=\"$DATASET_PREFIX\"|" | \
        sed "s|sample_filter=.*|sample_filter=\"$SAMPLE_FILTER\"|" | \
        sed "s|#SBATCH --job-name=.*|#SBATCH --job-name=CDSeq_${DATASET_PREFIX}_${SAMPLE_FILTER}|" > "$CDSEQ_SCRIPT"
    
    # Add custom bulk file parameter
    sed -i "s|Rscript.*|Rscript \${SCRIPT} \${RLIBRARY} \${input_data} \${output_base_dir} \${sample_filter} \"${BULK_FILE}\"|" "$CDSEQ_SCRIPT"
        
    chmod +x "$CDSEQ_SCRIPT"
    
    # Submit CDSeq job
    CDSEQ_JOB_ID=$(sbatch "$CDSEQ_SCRIPT" | grep -o '[0-9]*')
    if [ -n "$CDSEQ_JOB_ID" ]; then
        echo "Submitted CDSeq job $CDSEQ_JOB_ID" | tee -a "$MAIN_LOG"
        echo "CDSeq:${CDSEQ_JOB_ID}" >> "$GLOBAL_LOG_DIR/model_job_mapping.txt"
    else
        echo "ERROR: Failed to submit CDSeq job" | tee -a "$MAIN_LOG"
    fi
fi

echo "All deconvolution jobs submitted successfully." | tee -a "$MAIN_LOG"
echo "Results will be available in: $OUTPUT_DIR" | tee -a "$MAIN_LOG"