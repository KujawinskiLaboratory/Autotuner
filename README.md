# AutoTuner

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3270590.svg)](https://doi.org/10.5281/zenodo.3270590)
  [![Travis build status](https://travis-ci.org/crmclean/Autotuner.svg?branch=master)](https://travis-ci.org/crmclean/Autotuner)


## Introduction

This repo contains the code needed to run the R package AutoTuner. AutoTuner is used to identify dataset specific parameters to process untargeted metabolomics data. So far, AutoTuner has been tested on untargeted data generated on qTOF, orbitrap and Fourier transform ion cyclotron resonance mass analyzers. 

Currently, AutoTuner requires R version 3.4 or greater. 

For input, AutoTuner requires at least 3 samples of raw data converted from proprietary instrument formats (eg .mzML, .mzXML, or .CDF). It also requires a spreadsheet containing at least two columns. One column must match the raw data samples by name, and the other must describe the different experimental factors each sample belongs to. 

Please see vignettes/intro.Rmd for a tutorial on how to use AutoTuner within R. 

## AutoTuner Installation

The easiest way to use the package at the moment, is to download it using devtools. This can be accomplished by running the following code. 

```r
library(devtools)
install_github("crmclean/autotuner")
```

## mmetspData Installation 

Autotuner uses a second data package to test and demonstrate its functions. To download the package, please run the following within R.

```r
library(devtools)
install_github("crmclean/mmetspdata")
```

## Using AutoTuner

For a guide on how to use AutoTuner to find optimized data processing parameters, see vignettes/Autotuner.rmd
