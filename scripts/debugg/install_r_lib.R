#!/usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("2 arguments must be supplied (R_LIBRARY_PATH DO_INSTALL)",
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
  
  # Function to check and install packages with error handling
  install_if_missing <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cat(paste("Installing package:", pkg, "\n"))
      tryCatch({
        install.packages(pkg)
        if (requireNamespace(pkg, quietly = TRUE)) {
          cat(paste("Successfully installed package:", pkg, "\n"))
        } else {
          cat(paste("Failed to install package:", pkg, "\n"))
        }
      }, error = function(e) {
        cat(paste("Error installing package:", pkg, "-", e$message, "\n"))
      })
    } else {
      cat(paste("Package already installed:", pkg, "\n"))
    }
  }
  
  # Install devtools first (needed for GitHub packages)
  install_if_missing("devtools")
  
  # Install BiocManager (needed for Bioconductor packages)
  install_if_missing("BiocManager")
  
  # Install key Seurat dependencies first
  cat("Installing key Seurat dependencies...\n")
  seurat_dependencies <- c(
    "Rcpp",
    "rlang",
    "cowplot",
    "Matrix",
    "irlba",
    "RANN",
    "dplyr",
    "patchwork",
    "reticulate"
  )
  
  for (pkg in seurat_dependencies) {
    install_if_missing(pkg)
  }
  
  # Install required CRAN packages (alphabetically ordered)
  required_cran_packages <- c(
    "babelwhale",
    "biomaRt",
    "Biobase",
    "dplyr",
    "ggplot2",
    "gridExtra",
    "knitr",
    "Matrix",
    "Metrics",
    "readr",
    "reshape2",
    "RColorBrewer",
    "scales",
    "stringr",
    "tidyr",
    "tools"
  )
  
  # Install each package and print status
  for (pkg in required_cran_packages) {
    install_if_missing(pkg)
  }
  
  # Now install Seurat (which has many dependencies)
  cat("Installing Seurat...\n")
  install_if_missing("Seurat")
  
  # Test Singularity installation
  tryCatch({
    if (requireNamespace("babelwhale", quietly = TRUE)) {
      cat("Testing Singularity installation...\n")
      babelwhale::test_singularity_installation(detailed = TRUE)
    } else {
      cat("babelwhale not installed, skipping Singularity installation test\n")
    }
  }, error = function(e) {
    cat(paste("Error testing Singularity installation:", e$message, "\n"))
  })
  
  # Define a function to install from GitHub with error handling
  install_github <- function(repo, ...) {
    tryCatch({
      if (requireNamespace("devtools", quietly = TRUE)) {
        cat(paste("Installing", repo, "from GitHub...\n"))
        devtools::install_github(repo, ...)
        repo_name <- sub(".*/", "", repo)
        if (requireNamespace(repo_name, quietly = TRUE)) {
          cat(paste("Successfully installed", repo, "from GitHub\n"))
        } else {
          cat(paste("Failed to install", repo, "from GitHub\n"))
        }
      } else {
        cat("devtools not installed, skipping GitHub installation\n")
      }
    }, error = function(e) {
      cat(paste("Error installing", repo, "from GitHub:", e$message, "\n"))
    })
  }
  
  # Install GitHub packages
  cat("Installing GitHub packages...\n")
  install_github("tinnlab/DeconBenchmark")
  install_github("cartof/digitalDLSorter")
  install_github("xuranw/MuSiC")
  
  # Install CDSeq from CRAN or GitHub
  tryCatch({
    cat("Installing CDSeq from CRAN...\n")
    install.packages("CDSeq")
    if (requireNamespace("CDSeq", quietly = TRUE)) {
      cat("Successfully installed CDSeq from CRAN\n")
    } else {
      cat("Failed to install CDSeq from CRAN\n")
      # Try installing from GitHub if CRAN fails
      if (requireNamespace("devtools", quietly = TRUE)) {
        cat("Trying to install CDSeq from GitHub instead...\n")
        devtools::install_github("kkang7/CDSeq_R_Package")
      }
    }
  }, error = function(e) {
    cat(paste("Error installing CDSeq from CRAN:", e$message, "\n"))
    # Try installing from GitHub if CRAN fails
    if (requireNamespace("devtools", quietly = TRUE)) {
      cat("Trying to install CDSeq from GitHub instead...\n")
      tryCatch({
        devtools::install_github("kkang7/CDSeq_R_Package")
      }, error = function(e2) {
        cat(paste("Error installing CDSeq from GitHub:", e2$message, "\n"))
      })
    }
  })
  
  cat("Package installation completed.\n")
}