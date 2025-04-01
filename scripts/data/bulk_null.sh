#!/bin/sh
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 16G
#SBATCH --time 1:00:00
#SBATCH --output=/scratch/lorthiois/logs/%A.o
#SBATCH --error=/scratch/lorthiois/logs/%A.e
#SBATCH --job-name=bulk_null_gen

set -e  # Exit on error
set -x  # Print commands as they're executed

# Default parameters (can be overridden by command line arguments)
RLIBRARY="${1:-/work/gr-fe/R_4.3.1}"
SCRIPT="/work/gr-fe/lorthiois/DeconBenchmark/scripts/data/bulk_null.R"
INPUT_DATA="${2:-/work/gr-fe/lorthiois/DeconBenchmark/data/cleaned_feature_counts_matrix.csv}"
OUTPUT_BASE_DIR="${3:-/work/gr-fe/lorthiois/DeconBenchmark/generated_data}"
PREFIX="${4:-TB_D100}"

# Create prefix-specific subdirectory
OUTPUT_DIR="$OUTPUT_BASE_DIR/$PREFIX"
mkdir -p $OUTPUT_DIR

module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/
module load r

start=`date +%s`
echo "START AT $(date)"
echo "Generating null bulk data for prefix: $PREFIX"

# Make sure the R script exists
if [ ! -f "$SCRIPT" ]; then
  # If not found in the expected location, copy it there from the current directory
  if [ -f "bulk_null.R" ]; then
    echo "Copying bulk_null.R to $SCRIPT"
    cp bulk_null.R "$SCRIPT"
  else
    echo "ERROR: Could not find bulk_null.R script"
    exit 1
  fi
fi

# Run the R script
Rscript ${SCRIPT} ${RLIBRARY} ${INPUT_DATA} ${OUTPUT_DIR} ${PREFIX}

end=`date +%s`
runtime=$((end-start))
echo "Runtime: $runtime seconds"
echo "Null data saved to: $OUTPUT_DIR"