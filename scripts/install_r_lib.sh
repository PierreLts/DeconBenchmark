#!/bin/sh
## This script wil load R libraries
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 10G ##for 1
#SBATCH --time 6:00:00
#SBATCH --output=/scratch/lorthiois/logs/%A.o
#SBATCH --error=/scratch/lorthiois/logs/%A.e
#SBATCH --job-name=installRpackage


RLIBRARY="/work/gr-fe/R_4.3.1"
SCRIPT=/work/gr-fe/lorthiois/project2/scripts/install_r_lib.R
DOINSTALL=1 #if set to 1 will install package, if set to 0 will only check if it can load the libraries

module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/
module load r

eval "$(conda shell.bash hook)"
conda deactivate #Some package like RCurl will crah during installation if you don't deactive conda first

start=`date +%s`
echo "START AT $(date)"


Rscript ${SCRIPT} ${RLIBRARY} ${DOINSTALL}

# print end date and echo total runtime
end=`date +%s`
runtime=$((end-start))
echo Runtime: $runtime seconds
