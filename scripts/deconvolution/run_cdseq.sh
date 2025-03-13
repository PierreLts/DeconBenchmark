#!/bin/sh
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH --mem 64G
#SBATCH --time 12:00:00
#SBATCH --output=/scratch/lorthiois/logs/%A.o
#SBATCH --error=/scratch/lorthiois/logs/%A.e
#SBATCH --job-name=cdseq_run

### Useful for debugging 
set -e # Exit the script if any statement returns a non-true return value
set -x # Print each line of code being computed

RLIBRARY="/work/gr-fe/R_4.3.1" #IMPORTANT
SCRIPT=/work/gr-fe/lorthiois/DeconBenchmark/scripts/deconvolution/run_cdseq.R
dataset_prefix="${1:-TB}"  # Dataset prefix/subfolder
sample_filter="${2:-AB}"   # Sample filter: A, B, or AB (default: AB)
bulk_type="${3:-bulk}"     # New parameter: bulk file type (default: bulk)
input_data="/work/gr-fe/lorthiois/DeconBenchmark/generated_data/$dataset_prefix"  # Dataset directory
output_base_dir="/work/gr-fe/lorthiois/DeconBenchmark/deconv_results"

# Create dataset-specific output directory
output_dir="$output_base_dir/$dataset_prefix"
mkdir -p "$output_dir"

# Create and use temporary directory
TEMP_DIR="/scratch/lorthiois/data"
mkdir -p "$TEMP_DIR"
cd "$TEMP_DIR"
echo "Working in temporary directory: $TEMP_DIR"

module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/ #IMPORTANT
module load r #IMPORTANT

start=`date +%s`
echo "START AT $(date)"
echo "Processing dataset: $dataset_prefix with filter: $sample_filter"
echo "Using bulk file type: $bulk_type"

# Make sure CDSeq is installed
R --quiet --no-save << EOF
if (!requireNamespace("CDSeq", quietly = TRUE)) {
  install.packages("CDSeq", repos="https://cloud.r-project.org")
}
EOF

# Run the R script with the bulk_type parameter
Rscript ${SCRIPT} ${RLIBRARY} ${input_data} ${output_base_dir} ${sample_filter} ${bulk_type}

# Clean up temporary files
rm -f signature.csv bulk.csv
echo "Removed temporary files"

# Print end date and echo total runtime
end=`date +%s`
runtime=$((end-start))
echo "Runtime: $runtime seconds"
echo "Results saved to: $output_dir"