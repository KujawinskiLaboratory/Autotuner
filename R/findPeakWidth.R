#' @title findPeakWidth
#'
#' @description This function is designed to find the maximum peakwidth of an EIC observed within a
#' given TIC peak. It does so by using checkBounds to estimate width in time of a peak and countMaxima
#' to determine if a peak may be made up from two similar structural isomers.
#'
#' @param approvScorePeaks - A data.frame containing information on the
#' retained bins.
#' @param currentMsFile - mzR link to one of the entered raw data files.
#' @param sortedAllEIC - a da
#' @param boundaries - A numeric vector with indicies representing the scans
#' bounding the original TIC peak.
#' @param ppmEst - A scalar value representing the calculated ppm error
#' used to generate data.
#'
#' @return This function returns a scalar value representing an estimate for
#' the maximal peak width across samples.
#'
#' @export
findPeakWidth <- function(approvScorePeaks,
                          currentMsFile,
                          sortedAllEIC,
                          boundaries,
                          ppmEst) {

    maxScans <- max(approvScorePeaks$scanCount)
    maxPwTable <- approvScorePeaks[approvScorePeaks$scanCount == maxScans,]
    filteredRange <- unlist(maxPwTable[1,c("startMatch","endMatch")])

    ## finding the true bounds of the peak
    checkBoundaries <- filteredRange %in% boundaries

    sampleMetadata <- mzR::header(currentMsFile)
    sampleMetadata <- sampleMetadata[sampleMetadata$msLevel == 1L,]


    # Added this on 2019-03-24 for cases where ms2 data is not within the
    # ms convert file
    if(!all(sampleMetadata$msLevel == 1L)) {
        scans <- sub(".* scan=", "", sampleMetadata$spectrumId) %>% as.numeric()
    } else {
        scans <- 1:nrow(sampleMetadata)
    }

    ## adding this bandaid here to solve a problem I got with the FT data.
    ##
    if(max(filteredRange) > max(scans)) {
        scans <- sub(".* scan=", "", sampleMetadata$spectrumId) %>% as.numeric()
    }


    ## case 1 - there is a mz value spaning the range of the peak
    if(any(checkBoundaries)) {

        checkFeatures <- maxPwTable[order(maxPwTable$Intensity),]

        # Narrowing down number of features to 50 -----------------------------
        ## getting subset of masses to check for peak width
        if(nrow(checkFeatures) > 50) {
            checkMe <- round(sqrt(nrow(checkFeatures))) + 50

            if(checkMe > nrow(checkFeatures)) {
                checkMe <- nrow(checkFeatures)
            }

            checkTable <- checkFeatures[1:checkMe,]
            checkVals <- checkTable$meanMZ

        } else {
            checkTable <- checkFeatures
            checkVals <- checkTable$meanMZ
        }


        # checking the boundaries of peaks ------------------------------------
        ## looping through each of the features being checked
        print(length(checkVals))

        peakBounds <- list()
        for(massIndex in seq_along(checkVals)) {

            print(massIndex)
            mass <- checkVals[massIndex]
            # case 1 - the peak ends at both boundaries -----------------------
            if(length(checkBoundaries) == 2) {

                upperBound <- checkBounds(mass = mass,
                                          currentMsFile = currentMsFile,
                                          currentIndex = filteredRange[2],
                                          ppmEst = ppmEst,
                                          scans = scans,
                                          origBound = filteredRange[2],
                                          sampleMetadata = sampleMetadata)
                names(upperBound) <- "upper_bound"
                lowerBound <- checkBounds(mass = mass,
                                          upper = F,
                                          currentMsFile = currentMsFile,
                                          currentIndex = filteredRange[1],
                                          ppmEst = ppmEst,
                                          scans = scans,
                                          origBound = filteredRange[1],
                                          sampleMetadata = sampleMetadata)
                names(lowerBound) <- "lower_bound"
                peakBounds[[massIndex]] <- c(lowerBound, upperBound)

            } else {

            # the peak is observed a single boundary -------------------------


                # Case 1 - the peak is only bounded above ---------------------
                if(which(checkBoundaries) == 2) {

                    upperBound <- checkBounds(mass,
                                            currentMsFile = currentMsFile,
                                            currentIndex = filteredRange[2],
                                            ppmEst = ppmEst, sampleMetadata = sampleMetadata)
                    lowerBound <- checkTable$startMatch[massIndex]

                  ## case 2 - it is bounded from below
                } else {
                # Case 2 - the peak is only bounded below ---------------------

                    lowerBound <- checkBounds(mass,
                                            upper = F,
                                            currentMsFile = currentMsFile,
                                            currentIndex = filteredRange[1],
                                            ppmEst = ppmEst, sampleMetadata = sampleMetadata)
                    upperBound <- checkTable$endMatch[massIndex]

                }

                peakBounds[[massIndex]] <- c(lowerBound, upperBound)
            }

        }


        names(peakBounds) <- checkVals

        peakBounds <- unique(unlist(peakBounds))
        peakBounds <- peakBounds[!is.na(peakBounds)]

        if(length(peakBounds) == 1) {
            maxPw <- 0
        } else {
            rtUpper <- sampleMetadata$retentionTime[grep(paste0("scan=","\\b",
                                                                max(peakBounds), "\\b"),sampleMetadata$spectrumId)]
            rtLower <- sampleMetadata$retentionTime[grep(paste0("scan=","\\b",
                                                                min(peakBounds), "\\b"),sampleMetadata$spectrumId)]
            maxPw <- rtUpper - rtLower
        }


    ## case 2 - all peaks are bounded within the range of the calculated TIC peak.
    } else {

        curBounds <- unlist(maxPwTable[1,grep("Match", colnames(maxPwTable))])
        maxTime <- sampleMetadata$retentionTime[scans == max(curBounds)]
        minTime <- sampleMetadata$retentionTime[scans == min(curBounds)]
        maxPw <- maxTime - minTime

    }

    return(maxPw)

}
