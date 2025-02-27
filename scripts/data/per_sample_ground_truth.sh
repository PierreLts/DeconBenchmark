#!/bin/sh
## This script generates per-sample ground truth cell type proportions
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH --mem 64G
#SBATCH --time 1:00:00
#SBATCH --output=/scratch/lorthiois/logs/%A.o
#SBATCH --error=/scratch/lorthiois/logs/%A.e
#SBATCH --job-name=sample_ground_truth

set -e
set -x

RLIBRARY="/work/gr-fe/R_4.3.1" 
SCRIPT=/work/gr-fe/lorthiois/DeconBenchmark/scripts/data/per_sample_ground_truth.R
SEURAT_FILE="/work/gr-fe/lorthiois/DeconBenchmark/data/GFB-33245_HFKJMDSXC_2_scRNAseqWTATBseverityrun1_Seurat.rds"
OUTPUT_FILE="/work/gr-fe/lorthiois/DeconBenchmark/generated_data/per_sample_ground_truth_proportions.rda"

# Check if script exists
if [ ! -f "$SCRIPT" ]; then
  echo "ERROR: Script not found at $SCRIPT"
  exit 1
fi

# Check if Seurat file exists
if [ ! -f "$SEURAT_FILE" ]; then
  echo "ERROR: Seurat file not found at $SEURAT_FILE"
  exit 1
fi

# Check output directory
OUTPUT_DIR=$(dirname "$OUTPUT_FILE")
if [ ! -d "$OUTPUT_DIR" ]; then
  echo "Creating output directory: $OUTPUT_DIR"
  mkdir -p "$OUTPUT_DIR"
fi
if [ ! -w "$OUTPUT_DIR" ]; then
  echo "ERROR: Output directory is not writable: $OUTPUT_DIR"
  exit 1
fi

module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/
module load r

echo "R library path: $RLIBRARY"
echo "Script path: $SCRIPT"
echo "Seurat file: $SEURAT_FILE"
echo "Output file: $OUTPUT_FILE"

start=`date +%s`
echo "START AT $(date)"

Rscript ${SCRIPT} ${RLIBRARY} ${SEURAT_FILE} ${OUTPUT_FILE}
RSCRIPT_EXIT_CODE=$?
if [ $RSCRIPT_EXIT_CODE -ne 0 ]; then
  echo "ERROR: Rscript failed with exit code $RSCRIPT_EXIT_CODE"
  exit $RSCRIPT_EXIT_CODE
fi

end=`date +%s`
runtime=$((end-start))
echo "Runtime: $runtime seconds"