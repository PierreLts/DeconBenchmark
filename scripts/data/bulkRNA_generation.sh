#!/bin/sh
## This script wil launch a pQTL analysis
#the first time you run the script you'll need to set need_processing_genotype = TRUE, this will create the .bk and .rds file
#whenever you want to use the script on a different machine, you'll need to recreate these two files
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH --mem 16G ##for 1
#SBATCH --time 6:00:00
##SBATCH --mail-user=pierre.lorthiois@epfl.ch
##SBATCH --mail-type=END,FAIL 
#SBATCH --output=/scratch/lorthiois/logs/%A.o
#SBATCH --error=/scratch/lorthiois/logs/%A.e
#SBATCH --job-name=bulkRNA_generation

###Usefull for debuging 
set -e # Exit the script if any statement returns a non-true return value.
set -x #Print each line of code being computed
#set -u #Exit the script if a variable was not initialized


RLIBRARY="/work/gr-fe/R_4.3.1" #IMPORTANT
SCRIPT=/work/gr-fe/lorthiois/DeconBenchmark/scripts/data/bulkRNA_generation.R
input_data="/work/gr-fe/lorthiois/DeconBenchmark/data/cleaned_feature_counts_matrix.csv"
output_data="/work/gr-fe/lorthiois/DeconBenchmark/generated_data"

module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/ #IMPORTANT
module load r #IMPORTANT


start=`date +%s`
echo "START AT $(date)"

#(ORDER IS IMPORTANT)
Rscript ${SCRIPT} ${RLIBRARY} ${input_data} ${output_data}

# print end date and echo total runtime
end=`date +%s`
runtime=$((end-start))
echo Runtime: $runtime seconds