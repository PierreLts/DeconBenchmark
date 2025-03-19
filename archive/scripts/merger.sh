#!/bin/sh
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 32G
#SBATCH --time 0:30:00
#SBATCH --output=/scratch/lorthiois/logs/%A.o
#SBATCH --error=/scratch/lorthiois/logs/%A.e
#SBATCH --job-name=batch_merger

set -e
set -x

# Default parameters
RLIBRARY="/work/gr-fe/R_4.3.1"
SCRIPT="/work/gr-fe/lorthiois/DeconBenchmark/scripts/data/merger.R"
OUTPUT_PREFIX="TB"
FILTER_TYPE="AB"  # Default filter type (AB, A, or B)
BATCH1="TB1"
BATCH2="TB2"
BATCH3="TB3"
BATCH4="TB4"

# Override with command line arguments if provided
if [ "$1" != "" ]; then
    OUTPUT_PREFIX="$1"
fi
if [ "$2" != "" ]; then
    FILTER_TYPE="$2"
fi
# Additional batches could be specified but we're using the defaults above

# Validate filter type
if [ "$FILTER_TYPE" != "A" ] && [ "$FILTER_TYPE" != "B" ] && [ "$FILTER_TYPE" != "AB" ]; then
    echo "ERROR: Filter type must be 'A', 'B', or 'AB'. Got '$FILTER_TYPE'"
    exit 1
fi

module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/
module load r

start=`date +%s`
echo "START AT $(date)"
echo "Merging batches: $BATCH1 $BATCH2 $BATCH3 $BATCH4"
echo "Output prefix: $OUTPUT_PREFIX"
echo "Filter type: $FILTER_TYPE"

Rscript ${SCRIPT} ${RLIBRARY} ${OUTPUT_PREFIX} ${FILTER_TYPE} ${BATCH1} ${BATCH2} ${BATCH3} ${BATCH4}

end=`date +%s`
runtime=$((end-start))
echo "Runtime: $runtime seconds"
echo "Results saved to: /work/gr-fe/lorthiois/DeconBenchmark/generated_data/${OUTPUT_PREFIX}"