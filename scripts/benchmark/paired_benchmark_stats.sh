#!/bin/sh
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 32G
#SBATCH --time 00:15:00
#SBATCH --output=/scratch/lorthiois/logs/%A.o
#SBATCH --error=/scratch/lorthiois/logs/%A.e
#SBATCH --job-name=paired_benchmark

RLIBRARY="/work/gr-fe/R_4.3.1" 
SCRIPT=/work/gr-fe/lorthiois/DeconBenchmark/scripts/benchmark/paired_benchmark_stats.R
RESULTS_DIR="/work/gr-fe/lorthiois/DeconBenchmark/deconv_results"
PER_SAMPLE_GT="/work/gr-fe/lorthiois/DeconBenchmark/generated_data/per_sample_ground_truth_proportions.rda"
OUTPUT_DIR="/work/gr-fe/lorthiois/DeconBenchmark/benchmark_results/paired_benchmarks"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/
module load r

start=`date +%s`
echo "START AT $(date)"

Rscript ${SCRIPT} ${RLIBRARY} ${RESULTS_DIR} ${PER_SAMPLE_GT} ${OUTPUT_DIR}

end=`date +%s`
runtime=$((end-start))
echo "Runtime: $runtime seconds"