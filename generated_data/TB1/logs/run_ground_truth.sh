#!/bin/bash
#SBATCH --job-name=ground_truth
#SBATCH --output=/scratch/lorthiois/logs/%j.o
#SBATCH --error=/scratch/lorthiois/logs/%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1:00:00

module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/
module load r
Rscript /work/gr-fe/lorthiois/DeconBenchmark/scripts/data/ground_truth.R /work/gr-fe/R_4.3.1 /work/gr-fe/lorthiois/DeconBenchmark/generated_data/TB1 /work/gr-fe/lorthiois/DeconBenchmark/generated_data/TB1 TB1
