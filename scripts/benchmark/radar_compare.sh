#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=0:30:00
#SBATCH --output=/scratch/lorthiois/logs/%A.o
#SBATCH --error=/scratch/lorthiois/logs/%A.e
#SBATCH --job-name=radar_plot

# This script generates radar plots from benchmark files

#############################################################
# CONFIGURATION - MODIFY THESE LINES AS NEEDED
#############################################################
# List the benchmark files you want to include in the radar plot
PLOT_TITLE="Bulk vs. Pseudobulk vs. Random Bulk"

DATASETS_TO_PLOT=(
  "TB_D100-bulk_benchmark_AB_select-AB.csv"
  "TB_D100-pseudobulk_benchmark_AB_select-AB.csv"
  "TB_D100-bulk_random_benchmark_AB_select-AB.csv"
)

  # "TB_D100-bulk_benchmark_A_select-A.csv"
  # "TB_D100-bulk_benchmark_AB_select-A.csv"
  # "TB_D100-bulk_benchmark_B_select-B.csv"
  # "TB_D100-bulk_benchmark_AB_select-B.csv"

  # "TB_D100-bulk_benchmark_AB_select-AB.csv"
  # "TB_D100-pseudobulk_benchmark_AB_select-AB.csv"
  # "TB_D100-bulk_random_benchmark_AB_select-AB.csv"

    # "TB_D4-bulk_benchmark_AB_select-AB.csv"
    # "TB_D12-bulk_benchmark_AB_select-AB.csv"
    # "TB_D52-bulk_benchmark_AB_select-AB.csv"  
    # "TB_D100-bulk_benchmark_AB_select-AB.csv"
    # "TB_D300-bulk_benchmark_AB_select-AB.csv"
    # "TB_D500-bulk_benchmark_AB_select-AB.csv"
    # "TB_D1000-bulk_benchmark_AB_select-AB.csv"
    # "TB_D5000-bulk_benchmark_AB_select-AB.csv"
    # "TB_D10000-bulk_benchmark_AB_select-AB.csv"
    # "TB_D33348-bulk_benchmark_AB_select-AB.csv"


#############################################################

# Default parameters
RLIBRARY="/work/gr-fe/R_4.3.1"
OUTPUT_FILE="radar_comparison.pdf"
BENCHMARK_DIR="/work/gr-fe/lorthiois/DeconBenchmark/benchmark_results"

# Log start time
echo "Starting radar_compare.sh at $(date)"

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
      echo "Usage: $(basename $0) [options] [filename1] [filename2] ..."
      echo "Options:"
      echo "  -o, --output FILE      Output file name (default: radar_comparison.pdf)"
      echo "  -l, --library PATH     R library path (default: /work/gr-fe/R_4.3.1)"
      echo "  -d, --directory PATH   Benchmark directory (default: /work/gr-fe/lorthiois/DeconBenchmark/benchmark_results)"
      echo "  -t, --title TITLE      Plot title (default: 'Performance Metrics Comparison')"
      echo "  -h, --help             Show this help message"
      exit 0
      ;;
    *)
      break
      ;;
  esac
done

# Create benchmark_files array
benchmark_files=()

# Check if any benchmark files are provided via command line
if [ $# -eq 0 ]; then
  echo "No command-line arguments provided. Using predefined dataset list."
  # Use the predefined DATASETS_TO_PLOT array
  set -- "${DATASETS_TO_PLOT[@]}"
else
  echo "Using command-line arguments instead of predefined dataset list."
fi

# Make sure we have files to process
if [ $# -eq 0 ]; then
  echo "Error: No benchmark files provided"
  exit 1
fi

# Process each input file
for input_file in "$@"; do
  echo "Processing: $input_file"
  
  # Case 1: It's a full path that exists directly
  if [[ -f "$input_file" ]]; then
    benchmark_files+=("$input_file")
    echo "  Using direct file path: $input_file"
    continue
  fi
  
  # Case 2: Just search through all subdirectories in the benchmark dir
  found=false
  
  # Use find to locate the file regardless of directory structure
  while IFS= read -r found_file; do
    benchmark_files+=("$found_file")
    echo "  Found file: $found_file"
    found=true
  done < <(find "$BENCHMARK_DIR" -type f -name "$input_file" 2>/dev/null)
  
  # If the file wasn't found, try to get more specific
  if [ "$found" = false ]; then
    echo "  File not found directly. Attempting to extract dataset name..."
    
    # Extract approximate dataset prefix (TB_DXXX-something)
    dataset_prefix=$(echo "$input_file" | grep -oE '^[^_]+_[^_]+-[^_]+')
    
    if [ -n "$dataset_prefix" ]; then
      echo "  Looking in directories for dataset: $dataset_prefix"
      
      # Search more specifically in benchmark directories
      potential_dirs=(
        "$BENCHMARK_DIR/$dataset_prefix/benchmarks"
        "$BENCHMARK_DIR/$dataset_prefix"
      )
      
      for dir in "${potential_dirs[@]}"; do
        if [ -d "$dir" ]; then
          echo "  Checking directory: $dir"
          
          # Look for exact name match
          if [ -f "$dir/$input_file" ]; then
            benchmark_files+=("$dir/$input_file")
            echo "  Found exact match: $dir/$input_file"
            found=true
            break
          fi
          
          # Look for partial matches if exact match not found
          for f in "$dir"/*; do
            if [[ -f "$f" && "$f" == *"$input_file"* ]]; then
              benchmark_files+=("$f")
              echo "  Found partial match: $f"
              found=true
              break 2
            fi
          done
        fi
      done
    fi
  fi
  
  # If still not found, report it
  if [ "$found" = false ]; then
    echo "  WARNING: Could not find file matching '$input_file' in benchmark directories"
  fi
done

if [ ${#benchmark_files[@]} -eq 0 ]; then
  echo "Error: No valid benchmark files found"
  exit 1
fi

# Processing summary
echo "Processing ${#benchmark_files[@]} files for radar plot:"
for file in "${benchmark_files[@]}"; do
  echo "  - $file"
done

# Load R module
module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/
module load r

echo "Generating radar plot..."
# Run the R script
r_script="/work/gr-fe/lorthiois/DeconBenchmark/scripts/benchmark/radar_compare.R"

# Check if R script exists
if [ ! -f "$r_script" ]; then
  echo "Error: R script not found at $r_script"
  exit 1
fi

# Run R script
Rscript "$r_script" "$RLIBRARY" "$OUTPUT_FILE" "$PLOT_TITLE" "${benchmark_files[@]}"

# Check if the output file was created
if [ -f "$OUTPUT_FILE" ]; then
  echo "Radar plot created successfully: $OUTPUT_FILE"
else
  echo "Error: Failed to create radar plot"
  exit 1
fi