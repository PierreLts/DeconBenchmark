#!/bin/sh
## This script generates ground truth cell type proportions from single-cell labels
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 2
#SBATCH --mem 4G
#SBATCH --time 0:30:00
#SBATCH --output=/scratch/lorthiois/logs/%A.o
#SBATCH --error=/scratch/lorthiois/logs/%A.e
#SBATCH --job-name=ground_truth

set -e
set -x

RLIBRARY="/work/gr-fe/R_4.3.1" 
SCRIPT=/work/gr-fe/lorthiois/DeconBenchmark/scripts/data/ground_truth.R
INPUT_DATA="/work/gr-fe/lorthiois/DeconBenchmark/generated_data"
OUTPUT_FILE="/work/gr-fe/lorthiois/DeconBenchmark/generated_data/ground_truth_proportions.rda"

module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/
module load r

start=`date +%s`
echo "START AT $(date)"

Rscript ${SCRIPT} ${RLIBRARY} ${INPUT_DATA} ${OUTPUT_FILE}

end=`date +%s`
runtime=$((end-start))
echo "Runtime: $runtime seconds"