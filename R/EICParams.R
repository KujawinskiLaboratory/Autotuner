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
#'
#' @details The function CheckEICPeaks handles all the peak specific
#' computations.
#'
#' @export

EICparams <- function(Autotuner, massThresh, peak_table) {

    # Checking input ----------------------------------------------------------
    assertthat::assert_that(nrow(peak_table) > 0,
                          msg = "Peak table with 0 rows was entered into EICparams function.")

    # itterating between samples ----------------------------------------------
    totalEstimates <- list()
    for(j in unique(peak_table$Sample)) {

        message("Currently on sample ", j)
        currentTable <- peak_table[peak_table$Sample == j,]
        currentFile <- Autotuner@file_paths[j]

        # add API backend to make it work with netCDF files
        # this is relevant to mzR package
        if(grepl("CDF", currentFile)) {
            back <- "netCDF"
        } else {
            back <- "pwiz"
        }

        currentMsFile <- mzR::openMSfile(currentFile, backend = back)
        rm(currentFile)

        # going through each peak from a sample -----------------------------
        pickedParams <- list()
        for(curPeak in 1:nrow(currentTable)) {

            start <- currentTable[curPeak,"Start_time"]
            end <- currentTable[curPeak,"End_time"]
            width <- currentTable$peak_width[curPeak]
            observedPeak <- list(start = start, end = end)

            ## currently here
            estimatedPeakParams <- checkEICPeaks(currentMsFile = currentMsFile,
                                                 observedPeak = observedPeak,
                                                 massThresh)

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

        mzR::close(currentMsFile)
    }

    totalEstimates <- Reduce(rbind, totalEstimates)

    return(totalEstimates)
}

