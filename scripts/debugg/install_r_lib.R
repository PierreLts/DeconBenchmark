#!/usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("2 arguments must be supplied ()",
   call.=FALSE)
}


####### Parameter of main function
path_Rlibrary <- args[1]
do_install <- args[2]


.libPaths(path_Rlibrary, FALSE)
print(.libPaths())
print(getRversion())


if(do_install == 1){
  # Set the CRAN mirror URL
  cran_mirror <- "https://cloud.r-project.org"
  # Set the mirror option
  options(repos = c(CRAN = cran_mirror))
  
  # Function to check and install packages
  install_if_missing <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
  }
  library(devtools)
  # Install required libraries
  # install_if_missing("ggplot2")
  # install_if_missing("reshape2")
  # install_if_missing("corrplot")
  # install_if_missing("pheatmap")
  # install_if_missing("Metrics")
  # install_if_missing("stats")
  # install_if_missing("scales")
    # Install required libraries
  required_packages <- c("dplyr", "tidyr", "readr", "stringr", "gridExtra")
  
  lapply(required_packages, install_if_missing)
  # Add to Singularity container definition or add this to run script
  if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
  devtools::install_github("cartof/digitalDLSorter")
  # # babelwhale
  # install_if_missing("babelwhale")
  # babelwhale::test_singularity_installation(detailed = TRUE)
  #library(CDSeq)
  # # DeconBenchmark
  # install_if_missing("devtools")
  # devtools::install_github("tinnlab/DeconBenchmark")
  #install.packages("CDSeq")
  # # install.packages("devtools")
  # devtools::install_github("kkang7/CDSeq_R_Package")
  # # install.packages("devtools")
  # devtools::install_github("kkang7/CDSeq_R_Package", build_vignettes = TRUE)


}






  #install.packages("MatrixEQTL", lib=path_Rlibrary)
  #install.packages("bigsnpr", lib=path_Rlibrary)
  
 
 ##Bioconductor packages
  #if (!require("BiocManager", quietly = TRUE))
  #  install.packages("BiocManager")

  #BiocManager::install("GSVA", force = TRUE)
#BiocManager::install("MungeSumstats")
#BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh38#BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38")


#BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh37")
#BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")


# #necessary for install_github
#     install.packages("devtools", lib=path_Rlibrary)
#     library(devtools)
#     #Github packages
#     install_github("baolinwu/KATSP")
#     install_github('gastonstat/AssotesteR')
#     remotes::install_github("YuLab-SMU/clusterProfiler")

# #Download tar.gz file at https://cran.r-project.org/src/contrib/Archive/RVtests/
#     #install.packages("C:/Users/boutrys/OneDrive - UCL/DRAFT_ARTICLES/Statistical_framework_ensemble/To_publish_Plos_computational_biology/Response_to_reviewers/SIMULATION/functions/library_to_install/RVtests_1.2.tar.gz", repos=NULL, type="source")
#     install.packages("${path_Rlibrary}/RVtests_1.2.tar.gz", repos=NULL, type="source", lib=path_Rlibrary)

