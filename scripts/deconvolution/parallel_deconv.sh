#!/bin/bash
# parallel_deconv.sh - Run multiple deconvolution methods in parallel

# Default parameters
DEFAULT_DATASET_PREFIX="TB1"
DEFAULT_METHODS="AdRoit,ARIC,AutoGeneS,BayCount,BayesCCE,BayesPrism,BayICE,BisqueMarker,BisqueRef,BseqSC,CellDistinguisher,CIBERSORT,CPM,DAISM,debCAM,Deblender,DeCompress,deconf,DeconICA,DeconPeaker,DeconRNASeq,deconvSeq,DecOT,DeMixT,DESeq2,digitalDLSorter,DSA,dtangle,DWLS,EMeth,EPIC,FARDEEP,ImmuCellAI,LinDeconSeq,Linseed,MCPcounter,MethylResolver,MIXTURE,MOMF,MuSic,MySort,NITUMID,PREDE,quanTIseq,ReFACTor,RNA-Sieve,scaden,SCDC,spatialDWLS,TOAST"

# "AdRoit,ARIC,AutoGeneS,BayCount,BayesCCE,BayesPrism,BayICE,BisqueMarker,BisqueRef,BseqSC,CellDistinguisher,CIBERSORT,CPM,DAISM,debCAM,Deblender,DeCompress,deconf,DeconICA,DeconPeaker,DeconRNASeq,deconvSeq,DecOT,DeMixT,DESeq2,digitalDLSorter,DSA,dtangle,DWLS,EMeth,EPIC,FARDEEP,ImmuCellAI,LinDeconSeq,Linseed,MCPcounter,MethylResolver,MIXTURE,MOMF,MuSic,MySort,NITUMID,PREDE,quanTIseq,ReFACTor,RNA-Sieve,scaden,SCDC,spatialDWLS,TOAST"

# Parse command line arguments
DATASET_PREFIX="${1:-$DEFAULT_DATASET_PREFIX}"
METHODS="${2:-$DEFAULT_METHODS}"

# Set paths
SCRIPT_DIR="/scratch/lorthiois/scripts"
mkdir -p "$SCRIPT_DIR"
TEMPLATE_DIR="/work/gr-fe/lorthiois/DeconBenchmark/scripts/deconvolution"
OUTPUT_BASE_DIR="/work/gr-fe/lorthiois/DeconBenchmark/deconv_results"
OUTPUT_DIR="$OUTPUT_BASE_DIR/$DATASET_PREFIX"
BENCHMARK_DIR="/work/gr-fe/lorthiois/DeconBenchmark/benchmark_results/$DATASET_PREFIX"
GROUND_TRUTH="/work/gr-fe/lorthiois/DeconBenchmark/generated_data/${DATASET_PREFIX}/${DATASET_PREFIX}_ground_truth_proportions.rda"
GLOBAL_LOG_DIR="/work/gr-fe/lorthiois/DeconBenchmark/logs/${DATASET_PREFIX}"
mkdir -p ${GLOBAL_LOG_DIR}

# Create output directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$BENCHMARK_DIR"

# Log file setup
LOG_DIR="$OUTPUT_DIR/logs"
mkdir -p "$LOG_DIR"
MAIN_LOG="$LOG_DIR/parallel_deconv_$(date +%Y%m%d_%H%M%S).log"

echo "==== Starting parallel deconvolution pipeline $(date) ====" | tee -a "$MAIN_LOG"
echo "Dataset prefix: $DATASET_PREFIX" | tee -a "$MAIN_LOG"
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
    echo "Submitting deconvolution job for model: $MODEL" | tee -a "$MAIN_LOG"
    
    # Create a temporary job script for this model
    TEMP_SCRIPT="$SCRIPT_DIR/temp_${MODEL}_${DATASET_PREFIX}_deconv.sh"
    
    # Copy the template and modify parameters
    cat "$TEMPLATE_DIR/deconv_run.sh" | \
        sed "s|dataset_prefix=.*|dataset_prefix=\"$DATASET_PREFIX\"|" | \
        sed "s|deconv_method=.*|deconv_method=\"$MODEL\"|" | \
        sed "s|#SBATCH --job-name=.*|#SBATCH --job-name=${MODEL}_${DATASET_PREFIX}|" | \
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
    echo "${MODEL}:${JOB_ID}" >> "$GLOBAL_LOG_DIR/model_job_mapping.txt"
done

# Special case for CDSeq which requires a separate execution
if [[ "$METHODS" == *"CDSeq"* ]] || [[ "$METHODS" == *"All"* ]]; then
    echo "Submitting special CDSeq job" | tee -a "$MAIN_LOG"
    CDSEQ_SCRIPT="$SCRIPT_DIR/temp_cdseq_${DATASET_PREFIX}_deconv.sh"
    
    # Create CDSeq script by modifying the template
    cat "$TEMPLATE_DIR/run_cdseq.sh" | \
        sed "s|dataset_prefix=.*|dataset_prefix=\"$DATASET_PREFIX\"|" | \
        sed "s|#SBATCH --job-name=.*|#SBATCH --job-name=CDSeq_${DATASET_PREFIX}|" > "$CDSEQ_SCRIPT"
        
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

# Store job dependencies for evaluation
JOB_DEPENDENCY=$(IFS=,; echo "${DECONV_JOB_IDS[*]}")

# Submit evaluation job with dependency
EVAL_SCRIPT="$SCRIPT_DIR/evaluate_${DATASET_PREFIX}_results.sh"
cat > "$EVAL_SCRIPT" << EOF
#!/bin/bash
#SBATCH --job-name=eval_${DATASET_PREFIX}
#SBATCH --output=$LOG_DIR/%j.out
#SBATCH --error=$LOG_DIR/%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1:00:00

# Run evaluation and benchmarking scripts for $DATASET_PREFIX
echo "Running evaluation and benchmarking for dataset: $DATASET_PREFIX"
echo "Results directory: $OUTPUT_DIR"
echo "Benchmark output directory: $BENCHMARK_DIR"

# Here you would add commands to run your evaluation scripts
# For example:
# Rscript /path/to/evaluation_script.R $OUTPUT_DIR $BENCHMARK_DIR $DATASET_PREFIX

echo "Evaluation completed"
EOF

chmod +x "$EVAL_SCRIPT"
echo "Submitting evaluation job dependent on all deconvolution jobs" | tee -a "$MAIN_LOG"
EVAL_JOB_ID=$(sbatch --dependency=afterany:$JOB_DEPENDENCY "$EVAL_SCRIPT" | grep -o '[0-9]*')
echo "Submitted evaluation job $EVAL_JOB_ID - will run after all deconvolution jobs complete" | tee -a "$MAIN_LOG"

echo "All deconvolution jobs submitted successfully. Monitoring job completion." | tee -a "$MAIN_LOG"
echo "Results will be available in: $OUTPUT_DIR" | tee -a "$MAIN_LOG"