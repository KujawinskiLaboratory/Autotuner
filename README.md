# AutoTuner

<!-- badges: start -->
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3270590.svg)](https://doi.org/10.5281/zenodo.3270590)
  [![Travis build status](https://travis-ci.org/crmclean/Autotuner.svg?branch=master)](https://travis-ci.org/crmclean/Autotuner)
  [![Codecov test coverage](https://codecov.io/gh/crmclean/Autotuner/branch/master/graph/badge.svg)](https://codecov.io/gh/crmclean/Autotuner?branch=master)
  <!-- badges: end -->

## Our Paper is Out!

If you use AutoTuner, please cite us at:

McLean C. & Kujawinski L. B. (2020). "AutoTuner: High Fidelity and Robust Parameter Selection for Metabolomics Data Processing." _Analytical Chemistry_ https://doi.org/10.1021/acs.analchem.9b04804 

## Introduction

This repo contains the code needed to run the R package AutoTuner. AutoTuner is used to identify dataset specific parameters to process untargeted metabolomics data. So far, AutoTuner has been tested on untargeted data generated on qTOF, orbitrap and Fourier transform ion cyclotron resonance mass analyzers. 

Currently, AutoTuner requires R version 3.6 or greater. 

For input, AutoTuner requires at least 3 samples of raw data converted from proprietary instrument formats (eg .mzML, .mzXML, or .CDF). It also requires a spreadsheet containing at least two columns. One column must match the raw data samples by name, and the other must describe the different experimental factors each sample belongs to. 

## AutoTuner Installation

AutoTuner is now available through [bioconductor](https://bioconductor.org/packages/devel/bioc/html/Autotuner.html). The current released version of the package may be installed by running the following code. Please use the devel version of the package to obtain the most up to date version of the code. 

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("Autotuner")
```

The development version of the package may also be downloaded using devtools. This can be accomplished by running the following code. 

```r
library(devtools)
install_github("crmclean/autotuner")
```

## Using AutoTuner

For a guide on how to use AutoTuner to find optimized data processing parameters, see vignettes/Autotuner.rmd
