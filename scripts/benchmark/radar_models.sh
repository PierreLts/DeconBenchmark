#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=0:30:00
#SBATCH --output=/scratch/lorthiois/logs/%A.o
#SBATCH --error=/scratch/lorthiois/logs/%A.e
#SBATCH --job-name=radar_model_plot

# This script generates a radar plot showing all models from a single benchmark file

#############################################################
# CONFIGURATION - MODIFY THIS LINE AS NEEDED
#############################################################
# Specify the benchmark file to plot
DATASET_TO_PLOT="TB_D100-bulk_random_benchmark_AB_select-AB.csv"
#############################################################

# Default parameters
RLIBRARY="/work/gr-fe/R_4.3.1"
BENCHMARK_DIR="/work/gr-fe/lorthiois/DeconBenchmark/benchmark_results"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    -o|--output)
      OUTPUT_FILE="$2"
      shift 2
      ;;
    -l|--library)
      RLIBRARY="$2"
      shift 2
      ;;
    -d|--directory)
      BENCHMARK_DIR="$2"
      shift 2
      ;;
    -t|--title)
      PLOT_TITLE="$2"
      shift 2
      ;;
    -h|--help)
      echo "Usage: $(basename $0) [options] <benchmark_file>"
      echo "Options:"
      echo "  -o, --output FILE      Output file name (default: radar_models.pdf)"
      echo "  -l, --library PATH     R library path (default: /work/gr-fe/R_4.3.1)"
      echo "  -d, --directory PATH   Benchmark directory (default: /work/gr-fe/lorthiois/DeconBenchmark/benchmark_results)"
      echo "  -t, --title TITLE      Plot title (default: 'Model Performance Comparison')"
      echo "  -h, --help             Show this help message"
      exit 0
      ;;
    *)
      BENCHMARK_FILE="$1"
      shift
      ;;
  esac
done

# Check if benchmark file is provided via command line or use default
if [ -z "$BENCHMARK_FILE" ]; then
  echo "No benchmark file provided via command line. Using configured DATASET_TO_PLOT: $DATASET_TO_PLOT"
  BENCHMARK_FILE="$DATASET_TO_PLOT"
fi

# Extract the base name without "_benchmark" for plot title and output filename
BASE_NAME=$(basename "$BENCHMARK_FILE" .csv | sed 's/_benchmark//')
OUTPUT_FILE="radar_${BASE_NAME}.pdf"
PLOT_TITLE="Models for ${BASE_NAME}"

# Log start time
echo "Starting radar_models.sh at $(date)"
echo "Benchmark file: $BENCHMARK_FILE"
echo "Output file: $OUTPUT_FILE"
echo "Plot title: $PLOT_TITLE"

# Find the benchmark file if it's not a full path
if [ ! -f "$BENCHMARK_FILE" ]; then
  echo "File not found at direct path, searching in benchmark directory..."
  
  # Check if it's in the benchmark directory
  if [ -f "$BENCHMARK_DIR/$BENCHMARK_FILE" ]; then
    BENCHMARK_FILE="$BENCHMARK_DIR/$BENCHMARK_FILE"
    echo "Found at: $BENCHMARK_FILE"
  else
    # Try to search in subdirectories
    echo "Searching in subdirectories..."
    
    # Extract dataset prefix for better targeted search
    DATASET_PREFIX=$(echo "$BENCHMARK_FILE" | grep -oE '^[^_]+_[^_]+-[^_]+')
    
    if [ -n "$DATASET_PREFIX" ]; then
      echo "Looking in directories for dataset: $DATASET_PREFIX"
      
      # Search more specifically in benchmark directories
      POTENTIAL_DIRS=(
        "$BENCHMARK_DIR/$DATASET_PREFIX/benchmarks"
        "$BENCHMARK_DIR/$DATASET_PREFIX"
      )
      
      FOUND=false
      for DIR in "${POTENTIAL_DIRS[@]}"; do
        if [ -d "$DIR" ]; then
          echo "Checking directory: $DIR"
          
          # Look for exact name match
          if [ -f "$DIR/$BENCHMARK_FILE" ]; then
            BENCHMARK_FILE="$DIR/$BENCHMARK_FILE"
            FOUND=true
            echo "Found exact match: $BENCHMARK_FILE"
            break
          fi
          
          # Look for partial matches
          for F in "$DIR"/*; do
            if [[ -f "$F" && "$F" == *"$BENCHMARK_FILE"* ]]; then
              BENCHMARK_FILE="$F"
              FOUND=true
              echo "Found partial match: $BENCHMARK_FILE"
              break 2
            fi
          done
        fi
      done
      
      if [ "$FOUND" = false ]; then
        # Use find to search all subdirectories if specific search failed
        echo "Using find to search all subdirectories..."
        FOUND_FILE=$(find "$BENCHMARK_DIR" -type f -name "*$BENCHMARK_FILE*" -print -quit 2>/dev/null)
        
        if [ -n "$FOUND_FILE" ]; then
          BENCHMARK_FILE="$FOUND_FILE"
          echo "Found via find: $BENCHMARK_FILE"
        else
          echo "Error: Could not find benchmark file matching '$BENCHMARK_FILE'"
          exit 1
        fi
      fi
    else
      echo "Error: Could not extract dataset prefix from benchmark file name"
      exit 1
    fi
  fi
fi

# Check if file exists after all search attempts
if [ ! -f "$BENCHMARK_FILE" ]; then
  echo "Error: Benchmark file not found: $BENCHMARK_FILE"
  exit 1
fi

echo "Using benchmark file: $BENCHMARK_FILE"
echo "Output file: $OUTPUT_FILE"
echo "Plot title: $PLOT_TITLE"

# Load R module
module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/
module load r

echo "Generating radar plot..."
# Run the R script
R_SCRIPT="/work/gr-fe/lorthiois/DeconBenchmark/scripts/benchmark/radar_models.R"

# Check if R script exists
if [ ! -f "$R_SCRIPT" ]; then
  echo "Error: R script not found at $R_SCRIPT"
  exit 1
fi

# Run R script
Rscript "$R_SCRIPT" "$RLIBRARY" "$OUTPUT_FILE" "$PLOT_TITLE" "$BENCHMARK_FILE"

# Check if the output file was created
if [ -f "$OUTPUT_FILE" ]; then
  echo "Radar plot created successfully: $OUTPUT_FILE"
else
  echo "Error: Failed to create radar plot"
  exit 1
fi