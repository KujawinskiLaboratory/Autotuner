#' @title findPeakWidth
#'
#' @description This function is designed to find the maximum peakwidth of an EIC observed within a
#' given TIC peak. It does so by using checkBounds to estimate width in time of a peak and countMaxima
#' to determine if a peak may be made up from two similar structural isomers.
#'
#' @param approvScorePeaks - A data.frame containing information on the
#' retained bins.
#' @param mzDb - A list of data.frames containing the m/z and intensity values
#' from each scan's mass spectra.
#' @param header - A data.fame containing metadata on the sample like
#' spectra type (MS1 vs MS2), retention time, and scan count.
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
                          mzDb,
                          header,
                          sortedAllEIC,
                          boundaries,
                          ppmEst) {


    ## replacing traditional scan count index with actual scan index number
    ## 2019-06-19
    maxScans <- max(approvScorePeaks$scanCount)
    maxPwTable <- approvScorePeaks[approvScorePeaks$scanCount == maxScans,]
    filteredRange <- unlist(maxPwTable[1,c("startScan","endScan")])
    filteredRange[1] <- sortedAllEIC$scanID[filteredRange[1] == sortedAllEIC$scan][1]
    filteredRange[2] <- sortedAllEIC$scanID[filteredRange[2] == sortedAllEIC$scan][1]

    checkBoundaries <- filteredRange %in% boundaries
    #header <- header[header$msLevel == 1L,]

    # Added this on 2019-03-24 for cases where ms2 data is not within the
    # ms convert file
    # 2019-06-18 - ISSUE HERE REGARDING MATCHING SCANS TO IDS - ALSO RELATED TO DISSECT SCANS FUNCTION
    allScansInData <- as.numeric(sub(".* scan=", "", header$spectrumId))

    if(length(allScansInData) == 0) {
        stop("Error during findPeakWidth. allScansInData var is length 0. Check structure of raw data header file.")
    }

    ## CENSORED THIS 2019-06-19
    ## adding this bandaid here to solve a problem I got with the FT data.
    ##
    #if(max(filteredRange) > max(scans)) {
    #    scans <- sub("(.* )?scan=", "", header$spectrumId) %>% as.numeric()
    #}


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

        peakBounds <- list()
        for(massIndex in seq_along(checkVals)) {

            mass <- checkVals[massIndex]
            # case 1 - the peak ends at both boundaries -----------------------
            if(length(checkBoundaries) == 2) {

                upperBound <- checkBounds(mass = mass,
                                          mzDb = mzDb,
                                          currentIndex = filteredRange[2],
                                          ppmEst = ppmEst,
                                          scans = allScansInData,
                                          origBound = filteredRange[2],
                                          header = header)
                names(upperBound) <- "upper_bound"
                lowerBound <- checkBounds(mass = mass,
                                          upper = F,
                                          mzDb = mzDb,
                                          currentIndex = filteredRange[1],
                                          ppmEst = ppmEst,
                                          scans = allScansInData,
                                          origBound = filteredRange[1],
                                          header = header)
                names(lowerBound) <- "lower_bound"
                peakBounds[[massIndex]] <- c(lowerBound, upperBound)

            } else {

            # the peak is observed a single boundary -------------------------


                # Case 1 - the peak is only bounded above ---------------------
                if(which(checkBoundaries) == 2) {

                    upperBound <- checkBounds(mass,
                                            mzDb = mzDb,
                                            currentIndex = filteredRange[2],
                                            ppmEst = ppmEst,
                                            header = header,
                                            origBound = filteredRange[2],
                                            scans = allScansInData)
                    lowerBound <- checkTable$startMatch[massIndex]

                  ## case 2 - it is bounded from below
                } else {

                # Case 2 - the peak is only bounded below ---------------------

                    lowerBound <- checkBounds(mass,
                                            upper = F,
                                            mzDb = mzDb,
                                            currentIndex = filteredRange[1],
                                            ppmEst = ppmEst,
                                            header = header,
                                            origBound = filteredRange[1],
                                            scans = allScansInData)
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
            rtUpper <- header$retentionTime[grep(paste0("scan=","\\b",
                                                                max(peakBounds),
                                                                "\\b"),
                                                         header$spectrumId)]
            rtLower <- header$retentionTime[grep(paste0("scan=","\\b",
                                                                min(peakBounds),
                                                                "\\b"),
                                                         header$spectrumId)]
            maxPw <- rtUpper - rtLower
        }


    ## case 2 - all peaks are bounded within the range of the calculated TIC peak.
    } else {

        curBounds <- unlist(maxPwTable[1,grep("(start|end)Scan", colnames(maxPwTable))])
        maxTime <- header$retentionTime[scans == max(curBounds)]
        minTime <- header$retentionTime[scans == min(curBounds)]
        maxPw <- maxTime - minTime

    }

    return(maxPw)

}
