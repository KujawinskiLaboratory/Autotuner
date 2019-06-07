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
#' @param peak_table - table of peak width values extracted with the function
#' peak_width_table.
#' @param useGap - Parameter carried into checkEICPeaks that tells Autotuner
#' whether to use the gap statustic to determine the proper number of clusters
#' to use during ppm parameter estimation.
#' @param varExpThresh - Numeric value representing the variance explained
#' threshold to use if useGap is false.
#' @param returnPpmPlots - Boolean value that tells R to return plots for
#' ppm distributions.
#' @param plotDir - Path where to store plots.
#'
#' @details The function CheckEICPeaks handles all the peak specific
#' computations.
#'
#' @export

EICparams <- function(Autotuner, massThresh, peak_table, useGap = T,
                      varExpThresh = 0.8, returnPpmPlots = T,
                      plotDir = ".") {

    # Checking input ----------------------------------------------------------
    assertthat::assert_that(nrow(peak_table) > 0,
                          msg = "Peak table with 0 rows was entered into EICparams function.")

    if(returnPpmPlots) {
        if(!dir.exists(plotDir)) {
            message("Directory in plotDir did not exist. Using the current working directory instead.")
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
        msnObj <- suppressMessages(MSnbase::readMSData(files = currentFile, mode = "onDisk"))

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
            estimatedPeakParams <- checkEICPeaks(mzDb = mzDb,
                                                 header = header,
                                                 observedPeak = observedPeak,
                                                 massThresh,
                                                 useGap, varExpThresh,
                                                 returnPpmPlots, plotDir,
                                                 filename = basename(currentFile))

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

