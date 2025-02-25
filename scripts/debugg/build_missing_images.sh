#!/bin/bash
# save as build_missing_images.sh

# Create directory for local images
SING_IMG_DIR=/scratch/lorthiois/singularity_images
mkdir -p $SING_IMG_DIR

# List of images to build
IMAGES=(
  "deconvolution/bseqsc:latest"
  "deconvolution/celldistinguisher:latest"
  "deconvolution/cibersort:latest"
  "deconvolution/demixt:latest"
  "deconvolution/methylresolver:latest"
)

# Build images from Dockerfiles in DeconBenchmark-docker
for IMG in "${IMAGES[@]}"; do
  IMG_NAME=$(echo $IMG | cut -d/ -f2 | cut -d: -f1)
  echo "Building $IMG_NAME image..."
  singularity build --sandbox $SING_IMG_DIR/$IMG_NAME.sif docker-archive:///work/gr-fe/lorthiois/DeconBenchmark-docker/$IMG_NAME/Dockerfile
done