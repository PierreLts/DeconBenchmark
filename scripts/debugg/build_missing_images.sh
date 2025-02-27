#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=2:00:00
#SBATCH --job-name=fix_cache
#SBATCH --output=/scratch/lorthiois/logs/%A.o
#SBATCH --error=/scratch/lorthiois/logs/%A.e

### Usefull for debuging 
set -e # Exit the script if any statement returns a non-true return value.
set -x #Print each line of code being computed
#set -u #Exit the script if a variable was not initialized

# Create definition file directory
mkdir -p /scratch/lorthiois/def_files

# Create a definition file using debootstrap
cat > /scratch/lorthiois/def_files/test.def << EOF
Bootstrap: debootstrap
OSVersion: focal
MirrorURL: http://us.archive.ubuntu.com/ubuntu/

%post
    apt-get update
    apt-get install -y --no-install-recommends \
        r-base \
        r-base-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev
    
    # Install base R packages
    R --no-save -e "install.packages('Matrix', repos='https://cloud.r-project.org/')"
EOF

# Build the image
mkdir -p /scratch/lorthiois/singularity_images
singularity build /scratch/lorthiois/singularity_images/test.sif /scratch/lorthiois/def_files/test.def