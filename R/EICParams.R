#' @title EICparams
#'
#' @description This function is designed to calculate the recommended
#' parameters from EIC peaks. It is the main holder function for a lot of
#' different ones involved in calculating EIC parameters.
#'
#' @param Autotuner - An Autotuner objected containing sample specific raw
#' data.
#' @param massThresh - A generous exact mass error threshold used to estimate
#' PPM for features.
#' @param useGap - Parameter carried into checkEICPeaks that tells Autotuner
#' whether to use the gap statustic to determine the proper number of clusters
#' to use during ppm parameter estimation.
#' @param varExpThresh - Numeric value representing the variance explained
#' threshold to use if useGap is false.
#' @param returnPpmPlots - Boolean value that tells R to return plots for
#' ppm distributions.
#' @param plotDir - Path where to store plots.
#' @param verbose - Boolean value used to indicate whether checkEICPeaks
#' function returns messages to the console.
#'
#' @details The function CheckEICPeaks handles all the peak specific
#' computations.
#'
#' @return A data.frame of all peak specific estimates.
#'
#' @examples
#' library(devtools)
#' if(!require("mmetspData")) {
#'     install_github("crmclean/mmetspData")
#' }
#' library(mmetspData)
#' mmetspFiles <- c(system.file("mzMLs/mtab_mmetsp_ft_120815_24.mzML",
#' package = "mmetspData"), system.file("mzMLs/mtab_mmetsp_ft_120815_25.mzML",
#' package = "mmetspData"), system.file("mzMLs/mtab_mmetsp_ft_120815_26.mzML",
#' package = "mmetspData"))
#' metadata <- read.csv(system.file("mmetsp_metadata.csv",
#' package = "mmetspData"),stringsAsFactors = FALSE)
#' metadata <- metadata[metadata$File.Name %in%
#' sub(pattern = ".mzML", "",basename(mmetspFiles)),]
#' Autotuner <- createAutotuner(mmetspFiles, metadata,
#' file_col = "File.Name", factorCol = "Sample.Type")
#'
#' lag <- 25
#' threshold <- 3.1
#' influence <- 0.1
#'
#' signals <- lapply(getAutoIntensity(Autotuner),
#' ThresholdingAlgo, lag, threshold, influence)
#'
#' Autotuner <- isolatePeaks(Autotuner, returned_peaks = 10, signals)
#'
#' EICparams(Autotuner = Autotuner, massThresh = .005, verbose = FALSE,
#' returnPpmPlots = FALSE, useGap = TRUE)
#'
#' @export
EICparams <- function(Autotuner, massThresh, useGap = TRUE,
                      varExpThresh = 0.8, returnPpmPlots = TRUE,
                      plotDir = ".", verbose = TRUE) {

    peak_table <- Autotuner@peak_table

    # Checking input ----------------------------------------------------------
    assertthat::assert_that(nrow(peak_table) > 0,
                          msg = paste("Peak table with 0 rows was entered into",
                          "EICparams function."))

    if(returnPpmPlots) {
        if(!dir.exists(plotDir)) {
            message(paste("Directory in plotDir did not exist.",
                    "Using the current working directory instead."))
            plotDir <- "."
        }
    }


    # itterating between samples ----------------------------------------------
    totalEstimates <- list()
    for(j in unique(peak_table$Sample)) {

        message("Currently on sample ", j)
        currentTable <- peak_table[peak_table$Sample == j,]
        currentFile <- Autotuner@file_paths[j]

        # Adding msnbase functionality to replace mzR API
        msnObj <- suppressMessages(MSnbase::readMSData(files = currentFile,
                                                       mode = "onDisk",
                                                       msLevel. = 1))

        header <- suppressWarnings( MSnbase::header(msnObj))
        allMzs <- MSnbase::mz(msnObj)
        allInt <- MSnbase::intensity(msnObj)

        mzDb <- list()
        for(i in seq_along(allInt)) {
            mzDb[[i]] <- cbind(mz = allMzs[[i]],
                            intensity = allInt[[i]])
        }
        rm(allMzs, allInt, msnObj, i)

        # going through each peak from a sample -----------------------------
        pickedParams <- list()
        for(curPeak in 1:nrow(currentTable)) {

            message("--- Currently on peak: ", curPeak)
            start <- currentTable[curPeak,"Start_time"]
            end <- currentTable[curPeak,"End_time"]
            width <- currentTable$peak_width[curPeak]
            observedPeak <- list(start = start, end = end)

            ## currently here
            if(verbose) {
                estimatedPeakParams <- checkEICPeaks(mzDb = mzDb,
                                                     header = header,
                                                     observedPeak =
                                                         observedPeak,
                                                     massThresh,
                                                     useGap, varExpThresh,
                                                     returnPpmPlots, plotDir,
                                                     filename =
                                                         basename(currentFile))
            } else {
                estimatedPeakParams <- suppressMessages(checkEICPeaks(mzDb =
                                                                          mzDb,
                                                     header = header,
                                                     observedPeak =
                                                         observedPeak,
                                                     massThresh,
                                                     useGap, varExpThresh,
                                                     returnPpmPlots, plotDir,
                                                     filename =
                                                         basename(currentFile))
                                                     )
            }


            if(is.null(estimatedPeakParams)) {
                next
            }

            pickedParams[[curPeak]] <- cbind(estimatedPeakParams,
                                             startTime = start,
                                             endTime = end,
                                             sampleID = j)


        }

        sampleParams <- Reduce(rbind, pickedParams)

        totalEstimates[[j]] <- sampleParams

    }

    totalEstimates <- Reduce(rbind, totalEstimates)

    return(totalEstimates)
}

