#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --output=/scratch/lorthiois/logs/%A.o
#SBATCH --error=/scratch/lorthiois/logs/%A.e
#SBATCH --job-name=batch_patient_benchmarks

### Useful for debugging
set -e # Exit script on error
set -x # Print commands before execution

# Set paths
SCRIPT_DIR="/scratch/lorthiois/scripts"
mkdir -p "$SCRIPT_DIR"
TEMPLATE_DIR="/work/gr-fe/lorthiois/DeconBenchmark/scripts/benchmark"
RESULTS_DIR="/work/gr-fe/lorthiois/DeconBenchmark/deconv_results"
OUTPUT_DIR="/work/gr-fe/lorthiois/DeconBenchmark/benchmark_results/patient_specific"
GROUND_TRUTH="/work/gr-fe/lorthiois/DeconBenchmark/generated_data/ground_truth_proportions.rda"
PATIENT_GROUND_TRUTH_DIR="/work/gr-fe/lorthiois/DeconBenchmark/generated_data"
LOG_DIR="/work/gr-fe/lorthiois/DeconBenchmark/logs"

# Create output directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"

# Find all result files
result_files=$(find "$RESULTS_DIR" -name "results_*_*.rda")

# Process each result file
for result_file in $result_files; do
    # Extract method name
    method_name=$(basename "$result_file" | sed -E 's/results_([^_]+)_.+\.rda/\1/')
    
    # Create temporary job script
    temp_script="$SCRIPT_DIR/temp_${method_name}_patient_stats.sh"
    
    # Copy template and modify
    cat "$TEMPLATE_DIR/model_stats_patient.sh" | \
        sed "s|MEAN_GROUND_TRUTH=.*|MEAN_GROUND_TRUTH=\"$GROUND_TRUTH\"|" | \
        sed "s|PATIENT_GROUND_TRUTH_DIR=.*|PATIENT_GROUND_TRUTH_DIR=\"$PATIENT_GROUND_TRUTH_DIR\"|" | \
        sed "s|RESULTS=.*|RESULTS=\"$result_file\"|" | \
        sed "s|OUTPUT_DIR=.*|OUTPUT_DIR=\"$OUTPUT_DIR\"|" | \
        sed "s|#SBATCH --job-name=.*|#SBATCH --job-name=${method_name}_pstats|" > "$temp_script"
    
    chmod +x "$temp_script"
    
    # Submit job
    job_id=$(sbatch "$temp_script" | grep -o '[0-9]*')
    echo "Submitted job $job_id for patient-specific benchmarking of $method_name"
    
    # Record job mapping
    echo "${method_name}:${job_id}" >> "$LOG_DIR/patient_stats_mapping.txt"
done

echo "All patient-specific benchmarking jobs submitted"

# Submit comparison job to run after all individual benchmarks complete
COMPARE_SCRIPT="$SCRIPT_DIR/temp_compare_patient_models.sh"
cat > "$COMPARE_SCRIPT" << 'EOF'
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --output=/scratch/lorthiois/logs/%A.o
#SBATCH --error=/scratch/lorthiois/logs/%A.e
#SBATCH --job-name=compare_p_models

# Set paths
RLIBRARY="/work/gr-fe/R_4.3.1"
BENCHMARK_DIR="/work/gr-fe/lorthiois/DeconBenchmark/benchmark_results/patient_specific"
OUTPUT_DIR="/work/gr-fe/lorthiois/DeconBenchmark/benchmark_results/patient_specific/comparison"

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/
module load r

start=`date +%s`
echo "START AT $(date)"

# Run the R script to compare all method results
# TODO: Create a compare_patient_models.R script that reads all summary files
# and generates comparison visualizations

# Print end date and echo total runtime
end=`date +%s`
runtime=$((end-start))
echo "Runtime: $runtime seconds"
echo "Patient-specific model comparison results saved to: ${OUTPUT_DIR}"
EOF

chmod +x "$COMPARE_SCRIPT"
#sbatch "$COMPARE_SCRIPT"
echo "You can run the comparison script after all individual benchmarks complete"