#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --output=/scratch/lorthiois/logs/%A.o
#SBATCH --error=/scratch/lorthiois/logs/%A.e
#SBATCH --job-name=patient_stats

### Useful for debugging
set -e # Exit the script if any statement returns a non-true return value
set -x # Print each line of code being computed

RLIBRARY="/work/gr-fe/R_4.3.1" # R library path
SCRIPT="/work/gr-fe/lorthiois/DeconBenchmark/scripts/benchmark/model_stats_patient.R"
PATIENT_GROUND_TRUTH="/work/gr-fe/lorthiois/DeconBenchmark/generated_data/patient_ground_truth_proportions.rda"
RESULTS="/work/gr-fe/lorthiois/DeconBenchmark/deconv_results/results_BayesPrism_Batch1.rda" # Default, can be overridden
OUTPUT_DIR="/work/gr-fe/lorthiois/DeconBenchmark/benchmark_results/patient_specific"

# Extract method and data names from the results filename
RESULTS_FILENAME=$(basename "${RESULTS}")
METHOD=$(echo "${RESULTS_FILENAME}" | sed -E 's/results_([^_]+)_.+\.rda/\1/')
DATA=$(echo "${RESULTS_FILENAME}" | sed -E 's/results_[^_]+_(.+)\.rda/\1/')

# Create method-specific output directory
METHOD_OUTPUT_DIR="${OUTPUT_DIR}/${METHOD}_${DATA}"
mkdir -p ${METHOD_OUTPUT_DIR}

module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/
module load r

start=`date +%s`
echo "START AT $(date)"

# Run the R script (ORDER IS IMPORTANT)
Rscript ${SCRIPT} ${RLIBRARY} ${PATIENT_GROUND_TRUTH} ${RESULTS} ${METHOD_OUTPUT_DIR}

# Print end date and echo total runtime
end=`date +%s`
runtime=$((end-start))
echo "Runtime: $runtime seconds"
echo "Patient-specific benchmarking results saved to: ${METHOD_OUTPUT_DIR}"