# DeconvBenchmarking pipeline

## 1. Installations
To install all required libraries run 'sbatch install_r_lib.sh'

## 2. Data preparation
In the data folder, you hould have your .rds single cell data and .csv bulk data.
Modify with your correct data names and run 'sbatch scRNA_generation.sh', 'sbatch labelsRNA_generation.sh', and 'sbatch bulkRNA_generation.sh'

## 3. Deconvolution run
To deconvolute your bulk RNA- data, modify with your desired method and run 'sbatch deconv_run.sh'
Results will be saved under deconv_results

## 4. Benchmarking
run 'sbatch benchmarking.sh; with your specified results and ground truth data (scRNA-)

## 5. Next steps
1. Automate the benchmarking and data prep for several batches and methods.
2. Make a sorted list of performing methods.