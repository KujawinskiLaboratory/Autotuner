#' @title findTruePeaks
#'
#' @description This function is designed to filter out bins that don't come from
#' continuous scans. The idea is that after this stage, the data is ready for
#' parameter estimation.
#'
#' @param truePeaks - A lsit containing indicies representing each bin.
#' @param sortedAllEIC - All the raw ms1 data extracted from the EIC peak.
#'
#' @export
findTruePeaks <- function(truePeaks, sortedAllEIC) {


    # initializing storage ----------------------------------------------------
    ppmData <- list()
    counter <- 1


    # Checking if bin elements come from adj scans ----------------------------
    for(i in seq_along(truePeaks)) {

        pickedPeak <- truePeaks[[i]]

        peakData <- sortedAllEIC[pickedPeak,]
        peakData <- peakData[order(peakData$scanCounter),]

        # remove features that are duplicated.  -------------------------------
        if(nrow(peakData) == 1) {
            next()
        }

        ## checking to make sure features comes from adjacent scans
        scanDiff <- sort(unique(peakData$scanCounter)) %>% diff()

        ## checking that binned peaks:
        # are being picked up within consecutive scans
        peakDists <- (length(unique(scanDiff)) == 1 && unique(scanDiff) == 1)
        if(!peakDists) {
            next()
        }

        ## checking to see if any binned masses come from the same scan
        countsInScans <- table(peakData$scanCounter)
        moreInAScan <- any(as.vector(countsInScans) > 1)

        if(moreInAScan) {

            for(w in 1:(length(unique(peakData$scan)) - 1)) {

                curScan <- unique(peakData$scan)[w]
                nextScan <- unique(peakData$scan)[w+1]

                ### add condition here for whenever it is the end of the scan

                scanStates <- peakData[peakData$scan == curScan,]
                nextStates <- peakData[peakData$scan == nextScan,]
                peakData <- peakData[peakData$scan != nextScan & peakData$scan != curScan,]


                ## do this if there are two states in first scan
                if(w == 1) {
                    errorNext <- lapply(scanStates$mz, function(x) {
                        obsError <- abs(x - nextStates$mz)/x * 10^6
                    })

                    checkMins <- sapply(errorNext, which.min)

                    initialState <- sapply(seq_along(errorNext), function(x) {
                        errorNext[[x]][checkMins[x]]
                    }) %>% which.min()

                    scanStates <- scanStates[initialState,]
                }

                nextStateIndex <- sapply(scanStates$mz, function(x) {
                    obsError <- abs(x - nextStates$mz)/x * 10^6
                    intensityProb <- nextStates$intensity/sum(nextStates$intensity)
                    errorInverse <- 1/obsError
                    nextStateProb <- errorInverse/sum(errorInverse) * intensityProb
                    nextStateProb/sum(nextStateProb)
                }) %>% which.max()

                nextStates <- nextStates[nextStateIndex,]

                ## store new states
                bestStates <- rbind(scanStates,nextStates)
                peakData <- rbind(bestStates, peakData)

            }

            multipleInScan <- T

        } else {

            multipleInScan <- F
        }

        peakData <- peakData[order(peakData$scan),]

        obsPPM <- sapply(2:length(peakData$mz), function(mz) {
            estimatePPM(peakData$mz[(mz - 1)], peakData$mz[mz])
        })

        # storing output ----------------------------------------------------------
        ppmData[[counter]] <- data.frame(meanMZ = mean(peakData$mz),
                                         startScan = min(peakData$scan),
                                         endScan = max(peakData$scan),
                                         scanCount = length(peakData$scanCounter),
                                         Intensity = sum(peakData$intensity),
                                         meanIntensity = mean(peakData$intensity),
                                         intensityDispersion = sd(peakData$intensity),
                                         minIntensity = min(peakData$intensity),
                                         meanPPM = paste(signif(obsPPM), collapse = ";"),
                                         start = min(peakData$scanCounter),
                                         end = max(peakData$scanCounter),
                                         startMatch = min(peakData$dataMatchIndex),
                                         endMatch = max(peakData$dataMatchIndex),
                                         multipleInScan,
                                         stringsAsFactors = F)

        counter <- 1 + counter

    }
    approvedPeaks <- Reduce(rbind, ppmData)

    return(approvedPeaks)
}
