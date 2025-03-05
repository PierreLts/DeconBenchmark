#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=0:30:00
#SBATCH --output=/scratch/lorthiois/logs/%A.o
#SBATCH --error=/scratch/lorthiois/logs/%A.e
#SBATCH --job-name=model_report

# Define important paths
RLIBRARY="/work/gr-fe/R_4.3.1"
SCRIPT_PATH="/work/gr-fe/lorthiois/DeconBenchmark/scripts/debugg/generate_models_report.R"
OUTPUT_DIR="/work/gr-fe/lorthiois/DeconBenchmark/reports"

# Create output directory if needed
mkdir -p ${OUTPUT_DIR}

# Set Singularity temp directory
export SINGULARITY_TMPDIR=/scratch/lorthiois/temp
mkdir -p $SINGULARITY_TMPDIR
export SINGULARITYENV_APPEND_PATH=/scratch/lorthiois/temp

# Load R module
module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/
module load r

start=`date +%s`
echo "START AT $(date)"

# Run the R script
Rscript ${SCRIPT_PATH} ${RLIBRARY} ${OUTPUT_DIR}

# Print runtime
end=`date +%s`
runtime=$((end-start))
echo "Runtime: $runtime seconds"