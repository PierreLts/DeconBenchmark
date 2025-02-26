#!/bin/sh
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 32
#SBATCH --mem 128G ##for 1
#SBATCH --time 01:00:00
#SBATCH --output=/scratch/lorthiois/logs/%A.o
#SBATCH --error=/scratch/lorthiois/logs/%A.e
#SBATCH --job-name=build_images

### Usefull for debuging 
set -e # Exit the script if any statement returns a non-true return value.
set -x #Print each line of code being computed
#set -u #Exit the script if a variable was not initialized

# Create directory for local images
SING_IMG_DIR=/scratch/lorthiois/singularity_images
DOCKER_REPO_DIR=/work/gr-fe/R_4.3.1/DeconBenchmark-docker
TEMP_DIR=/scratch/lorthiois/temp_build
mkdir -p $SING_IMG_DIR
mkdir -p $TEMP_DIR

# Correct capitalization of directories
DIRECTORIES=(
  "BseqSC"
  "CellDistinguisher"
  "CIBERSORT" 
  "DeMixT"
  "MethylResolver"
)

# Lowercase image names
IMAGE_NAMES=(
  "bseqsc"
  "celldistinguisher"
  "cibersort"
  "demixt"
  "methylresolver"
)

echo "Starting image builds at $(date)"

# Build images using definition files
for i in "${!DIRECTORIES[@]}"; do
  DIR="${DIRECTORIES[$i]}"
  IMG="${IMAGE_NAMES[$i]}"
  
  echo "---------------------------------------------"
  echo "Building $IMG image from $DIR directory..."
  
  # Check if directory exists
  if [ ! -d "$DOCKER_REPO_DIR/$DIR" ]; then
    echo "ERROR: Directory not found at $DOCKER_REPO_DIR/$DIR"
    continue
  fi
  
  # Create a Singularity definition file from the Dockerfile
  DEF_FILE="$TEMP_DIR/Singularity.$IMG"
  echo "Creating Singularity definition file at $DEF_FILE"
  
  # Start with bootstrap from docker
  echo "Bootstrap: docker" > $DEF_FILE
  echo "From: rocker/r-ver:4.1.3" >> $DEF_FILE
  echo "" >> $DEF_FILE
  echo "%post" >> $DEF_FILE
  
  # If Dockerfile exists, extract installation commands
  if [ -f "$DOCKER_REPO_DIR/$DIR/Dockerfile" ]; then
    echo "Extracting installation commands from $DOCKER_REPO_DIR/$DIR/Dockerfile"
    
    # Extract RUN commands from Dockerfile and append to definition file
    grep -E "^RUN" "$DOCKER_REPO_DIR/$DIR/Dockerfile" | sed 's/^RUN /    /' >> $DEF_FILE
    
    # Also extract and adapt any COPY commands (as files might be needed)
    if grep -q "COPY" "$DOCKER_REPO_DIR/$DIR/Dockerfile"; then
      echo "" >> $DEF_FILE
      echo "# Files that would be copied in Docker" >> $DEF_FILE
      grep -E "^COPY" "$DOCKER_REPO_DIR/$DIR/Dockerfile" | sed 's/^COPY /# COPY: /' >> $DEF_FILE
    fi
  else
    echo "ERROR: Dockerfile not found at $DOCKER_REPO_DIR/$DIR/Dockerfile"
    echo "    echo 'No Dockerfile found, using basic R image'" >> $DEF_FILE
  fi
  
  # Add any additional setup needed for the container
  echo "" >> $DEF_FILE
  echo "%environment" >> $DEF_FILE
  echo "    export LC_ALL=C" >> $DEF_FILE
  
  echo "Building Singularity image from definition file..."
  singularity build $SING_IMG_DIR/${IMG}.sif $DEF_FILE
  
  # Verify if image was created
  if [ -f "$SING_IMG_DIR/${IMG}.sif" ]; then
    echo "Successfully built image: $SING_IMG_DIR/${IMG}.sif"
  else
    echo "Failed to build $IMG image"
  fi
  
  echo "Completed processing of $IMG"
  echo "---------------------------------------------"
done

echo "All build attempts complete at $(date)"
echo "Images available at $SING_IMG_DIR:"
ls -la $SING_IMG_DIR

# Clean up
rm -rf $TEMP_DIR