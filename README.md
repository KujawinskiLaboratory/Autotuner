# Autotuner

This repo contains the code needed to run the R package AutoTuner. AutoTuner is used to identify dataset specific parameters to process untargeted metabolomics data. So far, AutoTuner has been tested on untargeted data generated on qTOF, orbitrap and Fourier transform ion cyclotron resonance mass analyzers. 

Currently, AutoTuner requires R version 3.4 or greater. 

For input, AutoTuner requires at least 3 samples of raw data converted from proprietary instrument formats (eg .mzML, .mzXML, or .CDF). It also requires a spreadsheet containing at least two columns. One column must match the raw data samples by name, and the other must describe the different experimental factors each sample belongs to. 

Please see vignettes/intro.Rmd for a tutorial on how to use AutoTuner within R. 

