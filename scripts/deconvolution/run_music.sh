#!/bin/sh
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 32
#SBATCH --mem 128G
#SBATCH --time 12:00:00
#SBATCH --output=/scratch/lorthiois/logs/%A.o
#SBATCH --error=/scratch/lorthiois/logs/%A.e
#SBATCH --job-name=music_run

### Useful for debugging 
set -e # Exit the script if any statement returns a non-true return value
set -x # Print each line of code being computed

# Define default parameters
RLIBRARY="/work/gr-fe/R_4.3.1" #IMPORTANT
SCRIPT="/work/gr-fe/lorthiois/DeconBenchmark/scripts/deconvolution/run_music.R"
DEFAULT_DATASET_PREFIX="TB1"
dataset_prefix="${1:-$DEFAULT_DATASET_PREFIX}"  # Use default if not provided
input_data="/work/gr-fe/lorthiois/DeconBenchmark/generated_data/$dataset_prefix"  # Dataset directory
output_base_dir="/work/gr-fe/lorthiois/DeconBenchmark/deconv_results"

# Create dataset-specific output directory
output_dir="$output_base_dir/$dataset_prefix"
mkdir -p "$output_dir"

module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/ #IMPORTANT
module load r #IMPORTANT

start=`date +%s`
echo "START AT $(date)"
echo "Processing dataset: $dataset_prefix"

# Run the R script
Rscript ${SCRIPT} ${RLIBRARY} ${input_data} ${output_base_dir}

# Print end date and echo total runtime
end=`date +%s`
runtime=$((end-start))
echo "Runtime: $runtime seconds"
echo "Results saved to: $output_dir"