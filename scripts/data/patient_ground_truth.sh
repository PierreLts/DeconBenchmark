#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=/scratch/lorthiois/logs/%A.o
#SBATCH --error=/scratch/lorthiois/logs/%A.e
#SBATCH --job-name=patient_ground_truth

### Useful for debugging
set -e # Exit the script if any statement returns a non-true return value
set -x # Print each line of code being computed

# Set paths
RLIBRARY="/work/gr-fe/R_4.3.1"
SCRIPT_PATH="/work/gr-fe/lorthiois/DeconBenchmark/scripts/data/patient_ground_truth.R"
INPUT_DATA="/work/gr-fe/lorthiois/DeconBenchmark/data/GFB-33245_HFKJMDSXC_2_scRNAseqWTATBseverityrun1_Seurat.rds"
OUTPUT_DIR="/work/gr-fe/lorthiois/DeconBenchmark/generated_data"

# Make sure output directory exists
mkdir -p ${OUTPUT_DIR}

# Load modules
module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/
module load r

# Start timing
start=`date +%s`
echo "START AT $(date)"

# Run the R script
Rscript ${SCRIPT_PATH} ${RLIBRARY} ${INPUT_DATA} ${OUTPUT_DIR}

# Print end time and runtime
end=`date +%s`
runtime=$((end-start))
echo "Runtime: $runtime seconds"
echo "Patient-specific ground truth generation completed"