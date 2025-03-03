# DeconvBenchmarking pipeline
## 1. Installations
Installation of libraries and packages are done running 'install_r_lib.sh'

## 2. Data preparation
### Process bulk RNA-seq data
sbatch scripts/data/bulkRNA_generation.sh

### Process single-cell RNA-seq data
sbatch scripts/data/scRNA_generation.sh

### Process cell type labels
sbatch scripts/data/labelsRNA_generation.sh

### Combine all data into a single package
sbatch scripts/data/finalRDA_generation.sh

### Generate ground truth cell type proportions
sbatch scripts/data/ground_truth.sh

### Generate per-sample ground truth proportions
sbatch scripts/data/per_sample_ground_truth.sh


## 3. Deconvolution run
# Run all supported methods in parallel
sbatch scripts/deconvolution/parallel_deconv.sh [INPUT_RDA_FILE] [METHODS_COMMA_SEPARATED]

# Example with specific methods
sbatch scripts/deconvolution/parallel_deconv.sh data/Batch1.rda "MuSic,DWLS,CIBERSORT,BayesPrism"

# Special case for CDSeq (requires separate execution)
sbatch scripts/deconvolution/run_cdseq.sh



## 4. Visualization
# Generate plots comparing predicted vs. ground truth proportions
sbatch scripts/benchmark/per_sample_multi_plot.sh

# Generate performance metrics for all methods
sbatch scripts/benchmark/evaluate_results.sh

## 5. Next steps
