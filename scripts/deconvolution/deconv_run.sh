#!/bin/sh
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 64
#SBATCH --mem 256G ##for 1
#SBATCH --time 12:00:00
##SBATCH --mail-user=pierre.lorthiois@epfl.ch
##SBATCH --mail-type=END,FAIL 
#SBATCH --output=/scratch/lorthiois/logs/%A.o
#SBATCH --error=/scratch/lorthiois/logs/%A.e
#SBATCH --job-name=deconv_run

### Usefull for debuging 
set -e # Exit the script if any statement returns a non-true return value.
set -x #Print each line of code being computed
#set -u #Exit the script if a variable was not initialized


RLIBRARY="/work/gr-fe/R_4.3.1" #IMPORTANT
SCRIPT=/work/gr-fe/lorthiois/DeconBenchmark/scripts/deconvolution/deconv_run.R
input_data="/work/gr-fe/lorthiois/DeconBenchmark/data/Batch1.rda"
output_data="/work/gr-fe/lorthiois/DeconBenchmark/deconv_results"
deconv_method="MuSic"
sparse_conversion=FALSE # Enable sparse matrix conversion

module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/ #IMPORTANT
module load r #IMPORTANT


start=`date +%s`
echo "START AT $(date)"

export SINGULARITY_TMPDIR=/scratch/lorthiois/temp
mkdir -p $SINGULARITY_TMPDIR
export SINGULARITYENV_APPEND_PATH=/scratch/lorthiois/temp

#(ORDER IS IMPORTANT)
Rscript ${SCRIPT} ${RLIBRARY} ${input_data} ${output_data} ${deconv_method} ${sparse_conversion}

# print end date and echo total runtime
end=`date +%s`
runtime=$((end-start))
echo Runtime: $runtime seconds