#!/bin/sh
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH --mem 64G
#SBATCH --time 1:00:00
#SBATCH --output=/scratch/lorthiois/logs/%A.o
#SBATCH --error=/scratch/lorthiois/logs/%A.e
#SBATCH --job-name=merge_batches

set -e  # Exit on error
set -x  # Print commands as they're executed

# R library path
RLIBRARY="/work/gr-fe/R_4.3.1"

# Script path
SCRIPT="/work/gr-fe/lorthiois/DeconBenchmark/scripts/data/merge_batches.R"

# Output file path
OUTPUT_FILE="/work/gr-fe/lorthiois/DeconBenchmark/data/merged_batches.rds"

# Define paths for each batch (modify these according to your file paths)
BATCH1="/work/gr-fe/lorthiois/DeconBenchmark/data/GFB-33245_HFKJMDSXC_2_scRNAseqWTATBseverityrun1_Seurat.rds"
BATCH2="/work/gr-fe/lorthiois/DeconBenchmark/data/GFB-33251_HGVCGDSXC_1_scRNAseqWTATBseverityrun2_Seurat.rds"
BATCH3="/work/gr-fe/lorthiois/DeconBenchmark/data/GFB-34330_HGTTGDSXC_3_scRNAseqWTATBseverityrun3_Seurat.rds"
BATCH4="/work/gr-fe/lorthiois/DeconBenchmark/data/GFB-34331_HGTTGDSXC_4_scRNAseqWTATBseverityrun4_Seurat.rds"

# Create output directory if it doesn't exist
OUTPUT_DIR=$(dirname "$OUTPUT_FILE")
mkdir -p "$OUTPUT_DIR"

module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/
module load r

start=`date +%s`
echo "START AT $(date)"
echo "Merging Seurat objects in order (batch 1 on top of batch 2, etc.)"
echo "Batch 1: $BATCH1"
echo "Batch 2: $BATCH2" 
echo "Batch 3: $BATCH3"
echo "Batch 4: $BATCH4"

# Run the R script with the specified parameters
# Batches will be merged in the order they are listed here
Rscript ${SCRIPT} ${RLIBRARY} ${OUTPUT_FILE} ${BATCH1} ${BATCH2} ${BATCH3} ${BATCH4}

end=`date +%s`
runtime=$((end-start))
echo "Runtime: $runtime seconds"
echo "Merged Seurat object saved to: $OUTPUT_FILE"