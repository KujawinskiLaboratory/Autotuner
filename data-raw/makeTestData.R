library(testthat)
library(magrittr)
library(devtools)

if(!require("mmetspData")) {
    install_github("crmclean/mmetspData")
}
library(mmetspData)
library(Autotuner)

mmetspFiles <- c(system.file("mzMLs/mtab_mmetsp_ft_120815_24.mzML", package = "mmetspData"),
                 system.file("mzMLs/mtab_mmetsp_ft_120815_25.mzML", package = "mmetspData"),
                 system.file("mzMLs/mtab_mmetsp_ft_120815_26.mzML", package = "mmetspData"))

runfile <- read.csv(system.file("mmetsp_metadata.csv", package = "mmetspData"),
                    stringsAsFactors = F)

runfile <- runfile[runfile$File.Name %in% sub(pattern = ".mzML", "", basename(mmetspFiles)),]

## Loading Autotuner
Autotuner <- createAutotuner(mmetspFiles,
                             runfile,
                             file_col = "File.Name",
                             factorCol = "Sample.Type")

lag <- 20
threshold<- 3
influence <- 0.1
signal <- lapply(Autotuner@intensity, ThresholdingAlgo, lag, threshold, influence)


# setting table of peaks for each sample
peaks <- extract_peaks(Autotuner = Autotuner, returned_peaks = 10,
                       signals = signal)


peak_table <- peakwidth_table(Autotuner = Autotuner,
                              peakList = peaks,
                              returned_peaks = 7)
peak_difference <- peak_time_difference(peak_table)

