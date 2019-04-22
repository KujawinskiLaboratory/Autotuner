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

        ## added check to handle case where 2+ masses are observed in one scan
        masses <- peakData[order(peakData$scanCounter),"mz"]

        if(moreInAScan) {

            peakData$index <- 1:nrow(peakData)
            obsPPM <- list()
            for(w in peakData$index) {
                curScan <- peakData$scanCounter[w]
                neighborScans <- which(peakData$scanCounter == curScan + 1)

                adjPairs <- lapply(neighborScans, function(neighbor) {
                    c(peakData$mz[w], peakData$mz[neighbor])
                })

                obsPPM[[w]] <- sapply(adjPairs, function(pair) {
                    estimatePPM(pair[1], pair[2])
                }) %>% unlist()


            }
            obsPPM <- unlist(obsPPM)

            multipleInScan <- T

        } else {

            obsPPM <- sapply(2:length(masses), function(mz) {
                estimatePPM(masses[(mz - 1)], masses[mz])
            })

            multipleInScan <- F
        }


        # storing output ----------------------------------------------------------
        ppmData[[counter]] <- data.frame(meanMZ = mean(masses),
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
