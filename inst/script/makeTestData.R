library(devtools)
if(!require("mmetspData")) {
    install_github("crmclean/mmetspData")
}
library(mmetspData)
library(Autotuner)

mmetspFiles <- c(system.file("mzMLs/mtab_mmetsp_ft_120815_24.mzML",
                             package = "mmetspData"),
                 system.file("mzMLs/mtab_mmetsp_ft_120815_25.mzML",
                             package = "mmetspData"),
                 system.file("mzMLs/mtab_mmetsp_ft_120815_26.mzML",
                             package = "mmetspData"))

runfile <- read.csv(system.file("mmetsp_metadata.csv", package = "mmetspData"),
                    stringsAsFactors = FALSE)

runfile <- runfile[runfile$File.Name %in% sub(pattern = ".mzML", "",
                                              basename(mmetspFiles)),]

## Loading Autotuner
Autotuner <- createAutotuner(mmetspFiles,
                                runfile,
                                file_col = "File.Name",
                                factorCol = "Sample.Type")

#saveRDS(object = Autotuner, file = here::here("data/preSignalAuto.rds"))

lag <- 20
threshold<- 3
influence <- 0.1
signals <- lapply(getAutoIntensity(Autotuner),
                    ThresholdingAlgo, lag, threshold, influence)


Autotuner <- isolatePeaks(Autotuner, returned_peaks = 10, signals)


## object used to test whole dataset parameter return function
eicParamEsts <- EICparams(Autotuner = Autotuner,
                            massThresh = .005,
                            verbose = FALSE,
                            returnPpmPlots = FALSE,
                            useGap = TRUE)

#saveRDS(object = eicParamEsts,
#        file = here::here("inst/extdata/eicParamsEsts.rds"))
