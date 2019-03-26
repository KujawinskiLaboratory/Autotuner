# subFunctions ------------------------------------------------------------
#' @title estimatePPM
#'
#' @description This function provides a convenient way to measure ppm between
#' two exact masses.
#'
#' @param first - first mass entered
#' @param second - second mass entered
#'
#' @importFrom dplyr "%>%"
#' @importFrom grDevices "boxplot.stats"
#' @importFrom graphics "abline" "legend" "lines" "par" "plot"
#' @importFrom stats "acf" "cor" "dist" "filter" "kmeans" "lm" "median"
#'  "sd" "smooth.spline" "weighted.mean"
#' @importFrom utils "read.csv" "write.csv"
#'
#' @return The ppm error between the first and second entered mass
#'
#' @export
estimatePPM <- function(first, second) {
  abs(first-second)/first * 10 ^ 6
}

# Filtering Peaks from Noise ----------------------------------------------
#' @title filterPeaksfromNoise
#'
#' @description This function is desgined to perform the binning of potential
#' features. At this point in the algorithm, a potential feature is 2+ m/z
#' values that are within the generous exact mass error windown provided by
#' the user.
#'
#' @param matchedMasses - A data.frame containing information on the
#' retained bins.
#'
#' @return A list with entries for noise peaks and true peaks.
#'
#' @export
filterPeaksfromNoise <- function(matchedMasses) {

  ## Initializing variables for subsetting
  if(length(matchedMasses$values) %% 2 == 0) {
    list_length <- length(matchedMasses$values)/2
  } else {
    list_length <- as.integer(length(matchedMasses$values)/2) + 1
  }
  truePeaks <- vector("list", list_length)
  truePeakIndex <- 1
  lengthCounter <- 1
  lastTrue <- matchedMasses$values[1]

  ## subsets things that have a second mass w/in user error
  for(rleIndex in 1:length(matchedMasses$values)) {

    start <- lengthCounter

    if(lastTrue == T) {
      start <- start + 1
    }

    end <- lengthCounter + matchedMasses$lengths[rleIndex]  - 1

    if(isTRUE(matchedMasses$values[rleIndex])) {
      end <- end + 1
      lastTrue = T
      truePeaks[[truePeakIndex]] <- start:end
      truePeakIndex <- truePeakIndex + 1
    } else {
      lastTrue = F
      if(!exists("no_match")) {
        no_match <- start:end
      } else {
        no_match <- c(no_match, start:end)
      }

    }
    lengthCounter <- lengthCounter + matchedMasses$lengths[rleIndex]
  }

  if(is.null(truePeaks[[list_length]])) {
    truePeaks <- truePeaks[1:(list_length - 1)]
  }

  rm(list_length,
     truePeakIndex,
     lastTrue,
     lengthCounter,
     start,
     end,
     rleIndex)
  return(list(no_match, truePeaks))
}



# Grouping together Noise Data --------------------------------------------
#' @title estimateSNThresh
#'
#' @description This function is responsible for computing an estimated s/n
#' threshold.
#'
#' @param no_match - This is a vector of numerical indicies within the raw data
#' mapping to scan data considered to come from noise.
#' @param sortedAllEIC - This is the raw data from a single TIC peak.
#' @param approvedPeaks - This is a data.frame that contains information on
#' which peaks come from TIC data.
#'
#' @return returns an estimated s/n threshold value
#'
#' @export
estimateSNThresh <- function(no_match, sortedAllEIC, approvedPeaks) {

    noisePeakTable <- sortedAllEIC[no_match,]
    noise_noise <- noisePeakTable$intensity
    scanCount <- sortedAllEIC$scanCounter
    maxScan <- max(scanCount, na.rm = T)
    minScan <- min(scanCount, na.rm = T)
    scanCount <- scanCount[no_match]

    for(peakID in 1:nrow(approvedPeaks)) {

        peakStart <- approvedPeaks[peakID,"start"]
        peakEnd <- approvedPeaks[peakID,"end"]
        peakDist <- peakEnd - peakStart
        lowerBound <- peakStart - peakDist * 2

        if(lowerBound < minScan) {
            lowerBound <- minScan
        }

        upperBound <- peakEnd + peakDist * 2
        if(upperBound > maxScan) {
            upperBound <- maxScan
        }

        peakNoise <- noise_noise[which(scanCount >= lowerBound & scanCount <= upperBound)]
        fixedNoise <- peakNoise[-which(peakNoise %in% boxplot.stats(peakNoise)$out)]

        if(length(fixedNoise) == 0) {
            next()
        }

        fixedNoiseMean <- fixedNoise %>% mean(na.rm = T)
        fixedNoiseSD <- fixedNoise %>% sd(na.rm = T)
        Peak <- approvedPeaks$Intensity[peakID]

        if((Peak - fixedNoiseMean) > 3*fixedNoiseSD) {

            ## selects as true peak
            SigNoiseThresh <- (Peak - fixedNoiseMean)/fixedNoiseSD

            if(exists("SN")) {
                SN <- c(SN,SigNoiseThresh)
            } else {
                SN <- SigNoiseThresh
            }

        } else {
          next()
        }
    }

    if(!exists("SN")) {
        return(NA)
    }

    return(SN)
}


# Estimate ppm error ------------------------------------------------------
#' @title filterPpmError
#'
#' @description This function computes an estimate for the ppm error threshold.
#'
#' @param approvedPeaks - This is a data.frame with information on bins retained
#' after filtering with user input mz error threshold and continuity checks.
#' @param useGap - Parameter carried into checkEICPeaks that tells Autotuner
#' whether to use the gap statustic to determine the proper number of clusters
#' to use during ppm parameter estimation.
#' @param varExpThresh - Numeric value representing the variance explained
#' threshold to use if useGap is false.
#'
#' @details A distribution is created from the set of all ppm values identified.
#' The most dense peak of this distribution is assumed to represent the standard
#' ppm error of the data.
#'
#' @return This function returns a scalar value representing ppm error estimate.
#'
#' @export
filterPpmError <- function(approvedPeaks, useGap, varExpThresh) {

    ppmObs <- approvedPeaks$meanPPM
    ppmObs <- strsplit(split = ";", x = as.character(ppmObs)) %>%
        sapply(as.numeric) %>%
        unlist()
    totalPPM <- length(ppmObs)

    if(useGap) {
        gapStat <- cluster::clusGap(x = as.matrix(ppmObs),
                                    FUNcluster = kmeans,
                                    K.max = 5,
                                    B = 7,
                                    verbose = F)

        gapStat <- gapStat$Tab
        clusters <- min(which(diff(-gapStat[,3]) > 0)) + 1
        kmeansPPM <- kmeans(ppmObs, clusters)

    } else {

        ## estimating clustering based on hard coded 80% Vexp threshold
        clustCount <- 1
        varExp <- 0
        while(varExp < .8 && clustCount < length(ppmObs)/2) {
            kmeansPPM <- kmeans(ppmObs, clustCount)
            varExp <- kmeansPPM$betweenss/kmeansPPM$totss
            clustCount <- clustCount + 1
        }

    }

    ## cluster which contains smallest ppm values
    clusterSize <- table(kmeansPPM$cluster) %>% sort(decreasing = T)
    maxCluster <- names(clusterSize)[1]
    minCluster <- which(kmeansPPM$cluster == maxCluster)
    rm(clusterSize)

    x <- ppmObs
    n <- length(x)
    h <- 1
    gauss <- function(x) 1/sqrt(2*pi) * exp(-(x^2)/2)
    gaussDKE <- function(a, x) gauss((x - a)/h)/(n * h)

    bumps <- sapply(ppmObs[minCluster], gaussDKE, x)

    OutlierScore <- sapply(1:n, function(xx) {
        rowSums(bumps)[xx]/(sum(rowSums(bumps))/n)
    })

    scoreSub <- which(OutlierScore > 1)
    ppmEst <- max(ppmObs[scoreSub])
    ppmEst <- ppmEst + sd(ppmObs[scoreSub])*3

    return(ppmEst)
}

# Peak Width Estimation ---------------------------------------------------
#' @title countMaxima
#'
#' @description This function is designed to determine the maximum number of
#' scans of the longest peak in the sample. It does this using the signal
#' processing function used earlier. This function is used during the estimation
#' of the maximum peakwdith.
#'
#' @param peakBounds - What are the bounds of the peak.
#' @param sortedAllEIC - Original data object.
#' @param prevRange - Longest measured peak range.
#'
#' @return maxScans - the maximum number of scans across all peaks
#'
#' @export
countMaxima <- function(peakBounds, sortedAllEIC, prevRange) {

  maxScans <- list()
  for(peakIndex in 1:length(peakBounds)) {

    # Organizing data for comparison ------------------------------------------
    peak <- names(peakBounds)[peakIndex] %>% as.numeric()

    featInfoSub <- sortedAllEIC[estimatePPM(sortedAllEIC$mz, peak) < 1,]

    if(nrow(featInfoSub) == 0) {
      next
    }

    featInfoSub <- featInfoSub[order(featInfoSub$dataMatchIndex),]

    curRange <- prevRange[peakIndex,]

    curRange[,1] <= featInfoSub$dataMatchIndex && curRange[,2] >= featInfoSub$dataMatchIndex

    correctEntries <- sapply(featInfoSub$dataMatchIndex, function(entry) {
            curRange[,1] <= entry && curRange[,2] >= entry
          })

    ## This is useful for getting the correct set of intensity values for the identified peak
    featInfoSub <- featInfoSub[correctEntries,]

    # expading peak range to new identified interval --------------------------
    curBoundInfo <- peakBounds[[peakIndex]]
    intensity <- curBoundInfo[grep("intensity", names(curBoundInfo), ignore.case = T)]
    extendedBounds <- curBoundInfo[grep("index", names(curBoundInfo), ignore.case = T)] %>% unlist()

    totalIntensity <- c(intensity$lower_Intensity,
                        featInfoSub$intensity,
                        intensity$upper_Intensity)


    # Finding minima and lengths of peaks between minima ----------------------
    if(length(totalIntensity) >= 10) {

      lag <- (length(totalIntensity) / 5) %>% as.integer()
      threshold <- 1
      influence <- 1
      signalOutput <- ThresholdingAlgo(totalIntensity, lag, threshold, influence)
      signalPatterns <- rle(signalOutput$signals)

      maxScans[[peakIndex]] <- max(signalPatterns$lengths[signalPatterns$values != -1])

    } else {

      maxScans[[peakIndex]] <- diff(extendedBounds) %>% abs()

    }

  }

  if(length(maxScans) == 0) {
    maxScans[[1]] <- 0
  }

  return(maxScans)

}


