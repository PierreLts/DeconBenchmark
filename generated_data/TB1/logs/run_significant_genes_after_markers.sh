#!/bin/bash
#SBATCH --job-name=sig_genes_gen
#SBATCH --output=/scratch/lorthiois/logs/%j.o
#SBATCH --error=/scratch/lorthiois/logs/%j.e
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --dependency=afterok:Submitting job: markers_gen
Command: Rscript /work/gr-fe/lorthiois/DeconBenchmark/scripts/data/seurat_markers_generation.R /work/gr-fe/R_4.3.1 /work/gr-fe/lorthiois/DeconBenchmark/data/GFB-33245_HFKJMDSXC_2_scRNAseqWTATBseverityrun1_Seurat.rds /work/gr-fe/lorthiois/DeconBenchmark/generated_data/TB1 /work/gr-fe/lorthiois/DeconBenchmark/data/mart_export.txt TB1
Submitted job 32636666 for markers_gen
32636666

module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/
module load r
Rscript /work/gr-fe/lorthiois/DeconBenchmark/scripts/data/seurat_significant_genes_generation.R /work/gr-fe/R_4.3.1 "/work/gr-fe/lorthiois/DeconBenchmark/generated_data/TB1/TB1_markers.rda" /work/gr-fe/lorthiois/DeconBenchmark/generated_data/TB1 TB1
