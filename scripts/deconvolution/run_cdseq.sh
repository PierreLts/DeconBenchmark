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
input_data="/work/gr-fe/lorthiois/DeconBenchmark/data/Batch1.rda"
output_data="/work/gr-fe/lorthiois/DeconBenchmark/deconv_results"

module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/ #IMPORTANT
module load r #IMPORTANT

start=`date +%s`
echo "START AT $(date)"

# Make sure CDSeq is installed
R --quiet --no-save << EOF
if (!requireNamespace("CDSeq", quietly = TRUE)) {
  install.packages("CDSeq", repos="https://cloud.r-project.org")
}
EOF

# Create the script directory if it doesn't exist
mkdir -p $(dirname ${SCRIPT})

# Run the R script (ORDER IS IMPORTANT)
Rscript ${SCRIPT} ${RLIBRARY} ${input_data} ${output_data}

# Print end date and echo total runtime
end=`date +%s`
runtime=$((end-start))
echo "Runtime: $runtime seconds"