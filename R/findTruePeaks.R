#' @title findTruePeaks
#'
#' @description This function is designed to filter out bins that don't come
#' from continuous scans. The idea is that after this stage, the data is ready
#' for parameter estimation.
#'
#' @param truePeaks A list containing indicies representing each bin.
#' @param sortedAllEIC All the raw ms1 data extracted from the EIC peak.
#'
#' @return a list of candidate EIC regions
findTruePeaks <- function(truePeaks, sortedAllEIC) {
#

    # initializing storage ----------------------------------------------------
    ppmData <- list()
    counter <- 1

    # Checking if bin elements come from adj scans ----------------------------
    for(i in seq_along(truePeaks)) {

        pickedPeak <- truePeaks[[i]]

        peakData <- sortedAllEIC[pickedPeak,]
        peakData <- peakData[order(peakData$scan),]

        # remove features that are duplicated.  -------------------------------
        if(nrow(peakData) == 1) {
            next()
        }

        ## checking to make sure features comes from adjacent scans
        scanDiff <- sort(unique(peakData$scan))
        scanDiff <- diff(scanDiff)

        ## checking that binned peaks:
        # are being picked up within consecutive scans
        peakDists <- (length(unique(scanDiff)) == 1 && unique(scanDiff) == 1)
        if(!peakDists) {
            next()
        }

        ## checking to see if any binned masses come from the same scan
        countsInScans <- table(peakData$scan)
        moreInAScan <- any(as.vector(countsInScans) > 1)

        if(moreInAScan) {

            for(w in seq_len(length(unique(peakData$scan)) - 1)) {

                curScan <- unique(peakData$scan)[w]
                nextScan <- unique(peakData$scan)[w+1]

                ### add condition here for whenever it is the end of the scan

                scanStates <- peakData[peakData$scan == curScan,]
                nextStates <- peakData[peakData$scan == nextScan,]
                curObsRows <- peakData$scan == curScan | peakData$scan ==
                    nextScan

                peakData <- peakData[!curObsRows,]

                ## do this if there are two states in first scan

                if(w == 1) {

                    errorNext <- lapply(scanStates$mz, function(x) {
                        obsError <- abs(x - nextStates$mz)/x * 10^6
                    })

                    checkMins <- vapply(X = errorNext, FUN = which.min,
                                        FUN.VALUE = numeric(1))

                    initialState <- vapply(X = seq_along(errorNext),
                                           FUN = function(x) {
                        errorNext[[x]][checkMins[x]]
                    }, FUN.VALUE = numeric(1))
                    initialState <- which.min(initialState)

                    scanStates <- scanStates[initialState,]
                }

                nextStateIndex <- vapply(X = scanStates$mz,
                                         FUN = function(x) {

                    obsError <- abs(x - nextStates$mz)/x * 10^6


                    if(any(obsError == 0)) {

                        ## corner case - error is 0 and there are no other
                        ## options
                        if(length(x) == 1) {
                            obsError <- 0.001
                        } else {
                            obsError[obsError == 0] <-
                                min(obsError[obsError != 0])/10
                        }

                    }

                    intensityProb <- nextStates$intensity/
                        sum(nextStates$intensity)
                    errorInverse <- 1/obsError
                    nextStateProb <- errorInverse/sum(errorInverse) *
                        intensityProb
                    return(max(nextStateProb/sum(nextStateProb)))

                },
                FUN.VALUE = numeric(1))
                nextStateIndex <- which.max(nextStateIndex)

                nextStates <- nextStates[nextStateIndex,]

                ## store new states
                bestStates <- rbind(scanStates,nextStates)
                peakData <- rbind(bestStates, peakData)
                peakData <- peakData[order(peakData$scan),]

            }

            multipleInScan <- TRUE

        } else {

            multipleInScan <- FALSE
        }

        obsPPM <- vapply(X = 2:length(peakData$mz), FUN = function(mz) {
            estimatePPM(peakData$mz[(mz - 1)], peakData$mz[mz])
        }, FUN.VALUE = numeric(1))

        # storing output -------------------------------------------------------
        ppmData[[counter]] <- data.frame(meanMZ = mean(peakData$mz),
                                         startScan = min(peakData$scan),
                                         endScan = max(peakData$scan),
                                         scanCount = length(peakData$scan),
                                         Intensity = sum(peakData$intensity),
                                         meanIntensity = mean(
                                             peakData$intensity),
                                         intensityDispersion = sd(
                                             peakData$intensity),
                                         minIntensity = min(peakData$intensity),
                                         meanPPM = paste(signif(obsPPM),
                                                         collapse = ";"),
                                         multipleInScan,
                                         stringsAsFactors = FALSE,
                                         index = i)

        counter <- 1 + counter

    }
    approvedPeaks <- Reduce(rbind, ppmData)
    return(approvedPeaks)
}
