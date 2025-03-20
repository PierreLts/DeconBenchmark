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
# It will find benchmark files based on dataset names or use direct file paths

#############################################################
# CONFIGURATION - MODIFY THESE LINES TO CHANGE WHICH FILES TO PLOT
#############################################################
# List the benchmark files you want to include in the radar plot
# You can use:
#   - Full file paths: "/path/to/benchmark_file.csv"
#   - Dataset/file names: The script will look in $BENCHMARK_DIR/$DATASET/benchmarks/
#   - Dataset:selection: To find files ending with "_select-X.csv"
DATASETS_TO_PLOT=(
    # Examples - replace with your actual files:
    "TB_D100-bulk_benchmark_AB_select-A.csv"    # Blue line in your example  
    "TB_D100-bulk_benchmark_AB_select-B.csv"    # Pink line in your example
    "TB_D100-bulk_benchmark_AB_select-AB.csv"   # Green line in your example
    # Add more files here
)
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
    -h|--help)
      echo "Usage: $(basename $0) [options] <file1> <file2> ..."
      echo "Options:"
      echo "  -o, --output FILE      Output file name (default: radar_comparison.pdf)"
      echo "  -l, --library PATH     R library path (default: /work/gr-fe/R_4.3.1)"
      echo "  -d, --directory PATH   Benchmark directory (default: /work/gr-fe/lorthiois/DeconBenchmark/benchmark_results)"
      echo "  -h, --help             Show this help message"
      echo ""
      echo "You can provide full file paths or just dataset names."
      echo "If you provide dataset names, the script will look for files in benchmark directories."
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
  echo "No command-line arguments provided. Using predefined dataset list from script."
  # Use the predefined DATASETS_TO_PLOT array
  set -- "${DATASETS_TO_PLOT[@]}"
else
  echo "Using command-line arguments instead of predefined dataset list."
fi

# Make sure we have files to process
if [ $# -eq 0 ]; then
  echo "Error: No benchmark files provided in command line or in the script configuration"
  echo "Either provide files as arguments or modify the DATASETS_TO_PLOT array in the script"
  echo "Run $(basename $0) --help for usage information"
  exit 1
fi

for input_file in "$@"; do
  # Check if this is a direct file path
  if [[ -f "$input_file" ]]; then
    # It's a direct file path
    benchmark_files+=("$input_file")
    echo "Using file: $input_file"
    continue
  fi
  
  # Check for dataset:selection format
  if [[ "$input_file" == *":"* ]]; then
    dataset_name="${input_file%%:*}"
    selection="${input_file##*:}"
    
    # Look in the benchmark directory for any file ending with _select-X.csv
    benchmark_dir="${BENCHMARK_DIR}/${dataset_name}/benchmarks"
    
    if [ ! -d "$benchmark_dir" ]; then
      echo "Warning: Benchmark directory not found: $benchmark_dir"
      continue
    fi
    
    # Find files matching pattern
    matching_files=( "${benchmark_dir}"/*"_select-${selection}.csv" )
    
    if [ ${#matching_files[@]} -gt 0 ] && [ -f "${matching_files[0]}" ]; then
      benchmark_files+=("${matching_files[0]}")
      echo "Found: ${matching_files[0]}"
    else
      echo "Warning: No files found matching *_select-${selection}.csv in $benchmark_dir"
    fi
    continue
  fi
  
  # If it ends with .csv, look in the benchmark directory
  if [[ "$input_file" == *.csv ]]; then
    # Extract dataset from filename (assuming TB_X-Y format)
    dataset_base=$(echo "$input_file" | grep -oP '^[^_]+(_[^_]+)?')
    
    if [ -z "$dataset_base" ]; then
      echo "Warning: Could not extract dataset name from $input_file"
      continue
    fi
    
    # Check in both the benchmark directory and benchmarks subdirectory
    potential_paths=(
      "${BENCHMARK_DIR}/${dataset_base}/benchmarks/${input_file}"
      "${BENCHMARK_DIR}/${dataset_base}/${input_file}"
    )
    
    file_found=false
    for potential_path in "${potential_paths[@]}"; do
      if [ -f "$potential_path" ]; then
        benchmark_files+=("$potential_path")
        echo "Found: $potential_path"
        file_found=true
        break
      fi
    done
    
    if [ "$file_found" = false ]; then
      echo "Warning: Could not find $input_file in benchmark directories"
    fi
    continue
  fi
  
  # Otherwise, treat as simple dataset name and look for any benchmark file
  dataset_dir="${BENCHMARK_DIR}/${input_file}/benchmarks"
  if [ -d "$dataset_dir" ]; then
    # Find any benchmark CSV file
    benchmark_csv_files=( "${dataset_dir}"/*benchmark*.csv )
    
    if [ ${#benchmark_csv_files[@]} -gt 0 ] && [ -f "${benchmark_csv_files[0]}" ]; then
      benchmark_files+=("${benchmark_csv_files[0]}")
      echo "Found: ${benchmark_csv_files[0]} (first matching file)"
      echo "Note: Multiple files found, using the first one. To specify a particular file, use the full name."
    else
      echo "Warning: No benchmark files found for dataset $input_file"
    fi
  else
    echo "Warning: Dataset directory not found: $dataset_dir"
  fi
done

if [ ${#benchmark_files[@]} -eq 0 ]; then
  echo "Error: No valid benchmark files found"
  exit 1
fi

# Load R module
module use /work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core/
module load r

echo "Generating radar plot from ${#benchmark_files[@]} benchmark files..."
# Run R script with the same base name as this script
script_dir=$(dirname "$0")
r_script="${script_dir}/radar_compare.R"

# If R script doesn't exist at the same location, use a standard path
if [ ! -f "$r_script" ]; then
  r_script="/work/gr-fe/lorthiois/DeconBenchmark/scripts/benchmark/radar_compare.R"
fi

# Check if R script exists
if [ ! -f "$r_script" ]; then
  echo "Error: R script not found at $r_script"
  exit 1
fi

# Run R script
Rscript "$r_script" "$RLIBRARY" "$OUTPUT_FILE" "${benchmark_files[@]}"

# Check if the output file was created
if [ -f "$OUTPUT_FILE" ]; then
  echo "Radar plot created successfully: $OUTPUT_FILE"
else
  echo "Error: Failed to create radar plot"
  exit 1
fi