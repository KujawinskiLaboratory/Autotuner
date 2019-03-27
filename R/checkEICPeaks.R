#' @title checkEICPeaks
#'
#' @description This function is the outer most function used to check for
#' individual EIC peak specific parameters.
#'
#' @param currentMsFile - This is an mzR object containing sample specific
#' information to be loaded into Autotuner parameter estimation algorithm.
#' @param observedPeak - A list with names 'start' and 'end' containing
#' scalar values representing the calculated peak boundary points
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
#'
#' @export
checkEICPeaks <- function(currentMsFile,
                          observedPeak,
                          massThresh = 0.01,
                          useGap,
                          varExpThresh,
                          returnPpmPlots,
                          plotDir) {




    # extracting ms1 information for current peak -----------------------------
    suppressWarnings(rm(no_match,approvedPeaks))
    filename <- basename(currentMsFile@fileName)
    sampleChrom <- mzR::header(currentMsFile)

    scanDiff <- diff(sampleChrom$retentionTime[sampleChrom$msLevel == 1L])
    rate <- mean(scanDiff)
    rm(scanDiff)

    sortedAllEIC <- dissectScans(currentMsFile,
                                 observedPeak = observedPeak,
                                 sampleChrom = sampleChrom)


    assertthat::assert_that(nrow(sortedAllEIC) > 0,
                            msg = "Check dissectScans within checkEICPeaks. A table with 0 rows was returned. Make sure correct TIC peaks were selected.")
    boundaries <- range(sortedAllEIC$dataMatchIndex)


    # Checking if sorted scans pass mz error threshold ------------------------
    matchedMasses <- rle(diff(sortedAllEIC$mz) < massThresh)
    noiseAndPeaks <- filterPeaksfromNoise(matchedMasses)
    no_match <- noiseAndPeaks[[1]]
    truePeaks <- noiseAndPeaks[[2]]
    rm(noiseAndPeaks)

    approvedPeaks <- findTruePeaks(truePeaks, sortedAllEIC)
    message("-------- ", signif(nrow(approvedPeaks)/length(truePeaks)*100, digits = 2),
            " % of bins retained after findTruePeaks function.")
    message("-------- pre findTruePeaks bin count: " , length(truePeaks))
    message("-------- post findTruePeaks bin count: ", nrow(approvedPeaks))


    if(nrow(approvedPeaks) == 0) {
        message("No observed m/z value met was observed across adjacent scans below an error of massThresh.")
        return(NULL)
    }

    # Filtering data by variability and ppm checks ----------------------------
    ppmEst <- filterPpmError(approvedPeaks, useGap, varExpThresh,
                             returnPpmPlots, plotDir, observedPeak,
                             filename)
    rm(filename)
    assertthat::assert_that(!is.na(ppmEst),
                            msg = "Output of filterPpmError function was NA. Something may have gone wrong here with input.")

    ppmObs <- approvedPeaks$meanPPM
    ppmObs <- strsplit(split = ";", x = as.character(ppmObs)) %>%
        sapply(as.numeric)

    noisyBin <- lapply(ppmObs, function(ppm) {
        any(ppm > ppmEst)
    }) %>% unlist()
    approvScorePeaks <- approvedPeaks[!noisyBin,]

    # Estimating PeakPicking Parameters ---------------------------------------
    SNest <- estimateSNThresh(no_match,
                              sortedAllEIC, approvScorePeaks) %>% min()

    assertthat::assert_that(!is.na(SNest),
                            msg = "Output of estimateSNThresh within checkEICPeaks was NA. Something went wrong here.")

    scanEst <- approvScorePeaks$scanCount %>% min()

    ### Noise Intensity Estimate
    noiseEst <- approvScorePeaks$minIntensity %>%
        min() - 100

    if(noiseEst < 0) {
        noiseEst <- approvScorePeaks$minIntensity %>%
            min() + 10
    }

    ### Prefilter Intensity Estimate
    intensityEst <- approvScorePeaks$Intensity %>% min()/sqrt(2)
    if(intensityEst < 0) {
        intensityEst <- approvScorePeaks$Intensity %>% min()
    }

    ### peakWidth Estimate
    maxPw <- findPeakWidth(approvScorePeaks = approvScorePeaks,
                           currentMsFile = currentMsFile,
                           sortedAllEIC,
                           boundaries,
                           ppmEst)

    minPw <- scanEst * rate
    if(max(maxPw) < min(maxPw)) {

        message("One peak had a maximum peakwidth smaller than the minimum peakwidth estimate.")
        maxPw <- maxPw*2

    }


    # Returning Data ----------------------------------------------------------
    estimatedPeakData <- data.frame(ppm = ppmEst,
                                    noiseThreshold = noiseEst,
                                    peakCount = nrow(approvedPeaks),
                                    prefilterI = intensityEst,
                                    prefilterScan = scanEst,
                                    TenPercentQuanSN = unname(SNest),
                                    maxPw, minPw)

    return(estimatedPeakData)
}



