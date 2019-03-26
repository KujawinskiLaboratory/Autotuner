#' @title findTruePeaks
#'
#' @description This function is designed to filter out bins that either 1) have
#' two or more peaks originating from a single scans or 2) don't come from
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
    for(i in 1:length(truePeaks)) {

        pickedPeak <- truePeaks[[i]]

        peakData <- sortedAllEIC[pickedPeak,]
        peakData <- peakData[order(peakData$scanCounter),]

        # remove features that are duplicated.  -------------------------------
        if(nrow(peakData) == 1) {
            next()
        }

        ## filtering out peaks that are not made up of consecutive scans
        scanDiff <- sort(peakData$scanCounter) %>% diff()

        ## checking that binned peaks:
        # 1) dont come from the same sample
        # 2) are being picked up within consecutive scans
        peakDists <- (length(unique(scanDiff)) == 1 && unique(scanDiff) == 1)
        if(!peakDists) {
            next()
        }
        masses <- peakData[order(peakData$scanCounter),"mz"]

        obsPPM <- sapply(2:length(masses), function(mz) {
            estimatePPM(masses[(mz - 1)], masses[mz])
        })


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
                                         stringsAsFactors = F)

        counter <- 1 + counter

    }
    approvedPeaks <- Reduce(rbind, ppmData)

    return(approvedPeaks)
}
