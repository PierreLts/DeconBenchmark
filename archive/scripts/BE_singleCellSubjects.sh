#!/bin/sh
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 16G
#SBATCH --time 0:30:00
#SBATCH --output=/scratch/lorthiois/logs/%A.o
#SBATCH --error=/scratch/lorthiois/logs/%A.e
#SBATCH --job-name=BE_subjects

set -e  # Exit on error
set -x  # Print commands as they're executed

# Default parameters (can be overridden by command line arguments)
RLIBRARY="${1:-/work/gr-fe/R_4.3.1}"
OUTPUT_BASE_DIR="${2:-/work/gr-fe/lorthiois/DeconBenchmark/generated_data}"
PREFIX="${3:-BE}"  # Default prefix is BE, but can be overridden

# Create prefix-specific subdirectory
OUTPUT_DIR="$OUTPUT_BASE_DIR/$PREFIX"
mkdir -p $OUTPUT_DIR

# Path to R script
SCRIPT="/work/gr-fe/lorthiois/DeconBenchmark/scripts/data/BE_singleCellSubjects.R"

# Load R module
module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/
module load r

start=`date +%s`
echo "START AT $(date)"
echo "Using prefix: $PREFIX"

# Run R script
Rscript ${SCRIPT} ${RLIBRARY} ${OUTPUT_DIR} ${PREFIX}

end=`date +%s`
runtime=$((end-start))
echo "Runtime: $runtime seconds"
echo "Single cell subjects saved to: $OUTPUT_DIR with prefix: $PREFIX"