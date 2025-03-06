#!/bin/sh
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH --mem 64G
#SBATCH --time 6:00:00
#SBATCH --output=/scratch/lorthiois/logs/%A.o
#SBATCH --error=/scratch/lorthiois/logs/%A.e
#SBATCH --job-name=seurat_prep

set -e
set -x

RLIBRARY="/work/gr-fe/R_4.3.1"
SCRIPT=/work/gr-fe/lorthiois/DeconBenchmark/scripts/data/prepare_seurat_inputs.R
SEURAT_FILE="/work/gr-fe/lorthiois/DeconBenchmark/data/GFB-33245_HFKJMDSXC_2_scRNAseqWTATBseverityrun1_Seurat.rds"
OUTPUT_DIR="/work/gr-fe/lorthiois/DeconBenchmark/generated_data"
MAPPING_FILE="/work/gr-fe/lorthiois/DeconBenchmark/data/mart_export.txt"

module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/
module load r

start=`date +%s`
echo "START AT $(date)"
Rscript ${SCRIPT} ${RLIBRARY} ${SEURAT_FILE} ${OUTPUT_DIR} ${MAPPING_FILE}
end=`date +%s`
runtime=$((end-start))
echo Runtime: $runtime seconds