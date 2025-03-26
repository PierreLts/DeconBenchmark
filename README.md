# DeconvBenchmarking pipeline
This pipeline aims at running and benchmarking reference and reference-free deconvolution models with Tuberculosis RNA data.
The pipeline is build on the [tinnlab/DeconBenchmark](https://github.com/tinnlab/DeconBenchmark) package which unifies the running of 50 deconvolution models.

To this day, 19 models ran properly.

Pierre Lorthiois

## 1. Installations
Installation of libraries and packages are done running 'install_r_lib.sh'

## 2. Data preparation
Order of compilation is important for dependecies

Original .rds data found in the /data folder.
The compilable data is a .rda file composed of 2 matrixes (single cell and bulk) and one vector (labels).
The section formats our .rds file into the correct format.

*Process bulk RNA-seq data*
sbatch scripts/data/bulkRNA_generation.sh

*Process single-cell RNA-seq data*
sbatch scripts/data/scRNA_generation.sh
Isolates scRNA data in a matrix and maps gene names correctly using /data/mart_export.txt

*Process cell type labels*
sbatch scripts/data/labelsRNA_generation.sh

*Combine all data into a single package*
sbatch scripts/data/finalRDA_generation.sh

*Generate ground truth cell type proportions*
sbatch scripts/data/ground_truth.sh
Mean ground truth cell types proportion for quick performance evaluation

*Generate per-sample ground truth proportions*
sbatch scripts/data/per_sample_ground_truth.sh
Ground truth cell type proportion for benchmarking


## 3. Deconvolution run
*Run all supported methods in parallel*
sbatch scripts/deconvolution/parallel_deconv.sh [INPUT_RDA_FILE] [METHODS_COMMA_SEPARATED]

*Example with specific methods*
sbatch scripts/deconvolution/parallel_deconv.sh data/Batch1.rda "MuSic,DWLS,CIBERSORT,BayesPrism"

*Special case for CDSeq (requires separate execution)*
sbatch scripts/deconvolution/run_cdseq.sh



## 4. Visualization
*Generate plots comparing predicted vs. ground truth proportions*
sbatch scripts/benchmark/per_sample_multi_plot.sh

## 5. Debugging
*Generate error_summary.html file to ave a quick overview of errors durring deconvolution and plotting*
sbatch scripts/debugg/error_summary.sh
Extracts error patterns from log files found using job mapping files.

## 6. Archive
Results and script archive
