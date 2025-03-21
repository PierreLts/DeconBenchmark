#!/bin/sh
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 64
#SBATCH --mem 256G
#SBATCH --time 12:00:00
#SBATCH --output=/scratch/lorthiois/logs/%A.o
#SBATCH --error=/scratch/lorthiois/logs/%A.e
#SBATCH --job-name=deconv_run

### Useful for debugging 
set -e # Exit the script if any statement returns a non-true return value
set -x # Print each line of code being computed

RLIBRARY="/work/gr-fe/R_4.3.1" #IMPORTANT
SCRIPT=/work/gr-fe/lorthiois/DeconBenchmark/scripts/deconvolution/deconv_run.R
dataset_prefix="${1:-TB}"  # Dataset prefix/subfolder
sample_filter="${2:-AB}"   # Sample filter: A, B, or AB (default: AB)
output_base_dir="/work/gr-fe/lorthiois/DeconBenchmark/deconv_results"
deconv_method="${3:-MuSiC}"  # Default method, can be overridden
bulk_type="${4:-bulk}"      # Default bulk file type, can be overridden
sparse_conversion=FALSE # Enable sparse matrix conversion

# Create dataset-specific output directory
output_dir="$output_base_dir/$dataset_prefix"
mkdir -p "$output_dir"

module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/ #IMPORTANT
module load r #IMPORTANT

start=`date +%s`
echo "START AT $(date)"
echo "Processing dataset: $dataset_prefix"
echo "Sample filter: $sample_filter"
echo "Bulk file type: $bulk_type"
echo "Deconvolution method: $deconv_method"

export SINGULARITY_TMPDIR=/scratch/lorthiois/temp
mkdir -p $SINGULARITY_TMPDIR
export SINGULARITYENV_APPEND_PATH=/scratch/lorthiois/temp

#(ORDER IS IMPORTANT)
Rscript ${SCRIPT} ${RLIBRARY} ${dataset_prefix} ${sample_filter} ${output_base_dir}/${dataset_prefix}-${bulk_type} ${deconv_method} ${sparse_conversion} ${bulk_type}
# print end date and echo total runtime
end=`date +%s`
runtime=$((end-start))
echo "Runtime: $runtime seconds"
echo "Results saved to: $output_dir"