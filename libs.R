#!/usr/bin/env Rscript

if (!dir.exists("./lib")) {
  dir.create("./lib")
}

.libPaths("./lib")

install_if_missing <- function(packages) {
  for (package in packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
      install.packages(package, lib = "./lib", repos = "http://cran.us.r-project.org")
    }
  }
}

required_packages <- c("ggplot2", "pheatmap", "argparse", "grid", "gridExtra")
install_if_missing(required_packages)

install_bioconductor_if_missing <- function(packages) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", lib = "./lib", repos = "http://cran.us.r-project.org")
  }
  BiocManager::install(packages, lib = "./lib", update = FALSE)
}

bioconductor_packages <- c("WGCNA", "RCy3")
install_bioconductor_if_missing(bioconductor_packages)

library(WGCNA, lib.loc = "./lib")
library(ggplot2, lib.loc = "./lib")
library(pheatmap, lib.loc = "./lib")
library(RCy3, lib.loc = "./lib")
library(argparse, lib.loc = "./lib")
library(grid, lib.loc = "./lib")
library(gridExtra, lib.loc = "./lib")

cat("All required libraries are installed and loaded successfully.\n")