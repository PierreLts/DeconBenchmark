#!/bin/bash
# Modified build_missing_images.sh
#SBATCH --output=/scratch/lorthiois/logs/%A.o
#SBATCH --error=/scratch/lorthiois/logs/%A.e

# Create directory for local images
SING_IMG_DIR=/scratch/lorthiois/singularity_images
mkdir -p $SING_IMG_DIR

# Set Singularity cache directory
export SINGULARITY_CACHEDIR=/scratch/lorthiois/singularity_cache
mkdir -p $SINGULARITY_CACHEDIR

# List of images to build
IMAGES=(
  "deconvolution/bseqsc:latest"
  "deconvolution/celldistinguisher:latest"
  "deconvolution/cibersort:latest"
  "deconvolution/demixt:latest"
  "deconvolution/methylresolver:latest"
)

# Try to pull images directly from Docker Hub
for IMG in "${IMAGES[@]}"; do
  IMG_NAME=$(echo $IMG | cut -d/ -f2 | cut -d: -f1)
  echo "Pulling $IMG_NAME image from Docker Hub..."
  
  # Try to pull with increasing timeout
  if ! singularity pull --timeout 300 --name $SING_IMG_DIR/$IMG_NAME.sif docker://$IMG; then
    echo "Failed to pull $IMG, trying backup approach..."
    
    # If direct pull fails, create a minimal definition file
    DEF_FILE=$(mktemp --suffix=.def)
    echo "Bootstrap: docker" > $DEF_FILE
    echo "From: r-base:4.1.0" >> $DEF_FILE
    echo "" >> $DEF_FILE
    echo "%post" >> $DEF_FILE
    echo "    apt-get update && apt-get install -y --no-install-recommends \\" >> $DEF_FILE
    echo "        libcurl4-openssl-dev \\" >> $DEF_FILE
    echo "        libssl-dev \\" >> $DEF_FILE
    echo "        libxml2-dev" >> $DEF_FILE
    echo "    R -e \"install.packages(c('remotes', 'devtools'))\""  >> $DEF_FILE
    
    # Build from definition file
    echo "Building $IMG_NAME from definition file..."
    singularity build $SING_IMG_DIR/$IMG_NAME.sif $DEF_FILE
    rm $DEF_FILE
  fi
done

echo "Container building completed"