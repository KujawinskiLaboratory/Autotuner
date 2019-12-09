# subFunctions ------------------------------------------------------------
#' @title estimatePPM
#'
#' @description This function provides a convenient way to measure ppm between
#' two exact masses.
#'
#' @param first Numeric of length 1 respresenting the first mass entered
#' @param second Numeric of length 1 representing second mass entered
#'
#' @importFrom grDevices "boxplot.stats"
#' @importFrom graphics "abline" "legend" "lines" "par" "plot"
#' @importFrom stats "acf" "cor" "dist" "filter" "kmeans" "lm" "median"
#'  "sd" "smooth.spline" "weighted.mean"
#' @importFrom utils "read.csv" "write.csv"
#'
#' @return The ppm error between the first and second entered mass
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
#' @param matchedMasses A data.frame containing information on the
#' retained bins.
#'
#' @return A list with entries for noise peaks and true peaks.
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
    for(rleIndex in seq_along(matchedMasses$values)) {

        start <- lengthCounter

        if(lastTrue == TRUE) {
            start <- start + 1
        }

        end <- lengthCounter + matchedMasses$lengths[rleIndex]  - 1

        if(isTRUE(matchedMasses$values[rleIndex])) {
            end <- end + 1
            lastTrue = TRUE
            truePeaks[[truePeakIndex]] <- start:end
            truePeakIndex <- truePeakIndex + 1
        } else {
            lastTrue = FALSE
            if(!exists("no_match")) {
                no_match <- start:end
            } else {
                no_match <- c(no_match, start:end)
            }

        }
        lengthCounter <- lengthCounter + matchedMasses$lengths[rleIndex]
    }

    if(is.null(truePeaks[[list_length]])) {
        truePeaks <- truePeaks[seq_len(list_length - 1)]
    }

    rm(list_length,truePeakIndex,lastTrue,lengthCounter,start,end,rleIndex)
    return(list(no_match, truePeaks))
}



# Grouping together Noise Data --------------------------------------------
#' @title estimateSNThresh
#'
#' @description This function is responsible for computing an estimated s/n
#' threshold.
#'
#' @param no_match This is a vector of numerical indicies within the raw data
#' mapping to scan data considered to come from noise.
#' @param sortedAllEIC This is the raw data from a single TIC peak.
#' @param approvedPeaks This is a data.frame that contains information on
#' which peaks come from TIC data.
#'
#' @return returns an estimated s/n threshold value
estimateSNThresh <- function(no_match, sortedAllEIC, approvedPeaks) {

    noisePeakTable <- sortedAllEIC[no_match,]
    noise_noise <- noisePeakTable$intensity
    scanCount <- sortedAllEIC$scan
    maxScan <- max(scanCount, na.rm = TRUE)
    minScan <- min(scanCount, na.rm = TRUE)
    scanCount <- scanCount[no_match]

    ## generating index for subseting of fixed noise obj
    scanIntervals <- list()
    for(peakID in seq_len(nrow(approvedPeaks))) {

        peakStart <- approvedPeaks[peakID,"startScan"]
        peakEnd <- approvedPeaks[peakID,"endScan"]
        peakDist <- peakEnd - peakStart
        lowerBound <- peakStart - peakDist * 2

        if(lowerBound < minScan) {
            lowerBound <- minScan
        }

        upperBound <- peakEnd + peakDist * 2
        if(upperBound > maxScan) {
            upperBound <- maxScan
        }

        scanIntervals[[peakID]] <- which(minScan:maxScan %in% c(lowerBound,
                                                                upperBound))


    }

    ## calculating all fixed noise values
    scanRange <- (minScan:maxScan)
    fixedNoiseList <- list()
    counter <- 1
    for(scanID in seq_along(scanRange)) {

        peakNoise <- noise_noise[scanCount == scanRange[scanID]]

        if(length(peakNoise) == 0) {
            next
        }

        fixedNoise <- peakNoise[!(peakNoise %in% boxplot.stats(peakNoise)$out)]

        fixedNoiseMean <- mean(x = fixedNoise, na.rm = TRUE)
        fixedNoiseVar <-  stats::var(x = fixedNoise, na.rm = TRUE)

        ## 2019-07-08: fixed corner case where only one noise element was
        ## found within the bin
        if(is.na(fixedNoiseVar)) {
            fixedNoiseVar <- 0
        }

        N <- length(fixedNoise)

        fixedNoiseList[[counter]] <- data.frame(fixedNoiseMean,fixedNoiseVar,N)
        counter <- 1 + counter
    }
    fixedNoiseList <- Reduce(rbind, fixedNoiseList)

    ## calculating the sd and mean for each group
    noiseIntDb <- list()
    counter <- 1
    for(row in seq_along(scanIntervals)) {

        if(row == 1) {

            new <- TRUE

        } else {

            new <- vapply(X = noiseIntDb, FUN = function(noiseDb) {
                if(!all(scanIntervals[[row]] == as.numeric(noiseDb[,c(1,2)]))) {
                    return(TRUE)
                } else {
                    return(FALSE)
                }
            }, FUN.VALUE = logical(1))
            new <- all(new)
        }


        if(new) {
            ## add check to see if cur row has already been looked at
            curRow <- scanIntervals[[row]]
            curStatDb <- fixedNoiseList[curRow[1]:curRow[2],]

            ## added this to check for case when row is missing
            ## if a match to the noise peak was not identified
            if(any(is.na(curStatDb$fixedNoiseMean))) {
                curStatDb <- curStatDb[!is.na(curStatDb$fixedNoiseMean),]
            }

            eX2 <- sum((curStatDb$fixedNoiseMean^2 +
                            curStatDb$fixedNoiseVar)*curStatDb$N)
            eX2 <- eX2/sum(curStatDb$N)
            groupMean <- mean(curStatDb$fixedNoiseMean)
            groupVar <- eX2 - groupMean^2

            if(groupVar < 0) {
                groupVar <- groupVar * -1
            }

            groupSd <- suppressWarnings(sqrt(groupVar))
            noiseIntDb[[counter]] <- data.frame(start = curRow[1],
                                                end = curRow[2], groupMean,
                                                groupSd)
            counter <- counter + 1
        }

    }

    noiseIntDb <- Reduce(rbind, noiseIntDb)
    noiseIntDb$key <- apply(noiseIntDb[,c(1,2)], 1, paste, collapse = " ")
    rm(curStatDb, eX2, groupSd, groupVar, groupMean,
        curRow, fixedNoiseList)


    SN <- list()
    counter <- 1
    for(peakID in seq_along(scanIntervals)) {

        scanInt <- paste(scanIntervals[[peakID]], collapse = " ")
        scanStats <- noiseIntDb[noiseIntDb$key == scanInt,]

        if(nrow(scanStats) == 0) {
            next()
        }

        ## 2019-06-20 - added here to fix bug if noise calc is wrong
        if(is.nan(scanStats$groupSd) | is.na(scanStats$groupSd)) {
            next()
        }

        if(is.nan(scanStats$groupMean)) {
            next()
        }

        Peak <- approvedPeaks$Intensity[peakID]

        if((Peak - scanStats$groupMean) > 3*scanStats$groupSd) {

            ## selects as true peak
            SigNoiseThresh <- (Peak - scanStats$groupMean)/scanStats$groupSd
            SN[[counter]] <- SigNoiseThresh
            counter <- counter + 1

        } else {

            next()

        }
    }

    if(!exists("SN")) {
        return(NA)
    }

    return(unlist(SN))
}


# Estimate ppm error ------------------------------------------------------
#' @title filterPpmError
#'
#' @description This function computes an estimate for the ppm error threshold.
#'
#' @param approvedPeaks This is a data.frame with information on bins retained
#' after filtering with user input mz error threshold and continuity checks.
#' @param useGap Parameter carried into checkEICPeaks that tells Autotuner
#' whether to use the gap statustic to determine the proper number of clusters
#' to use during ppm parameter estimation.
#' @param varExpThresh Numeric value representing the variance explained
#' threshold to use if useGap is false.
#' @param returnPpmPlots Boolean value that tells R to return plots for
#' ppm distributions.
#' @param plotDir Path where to store plots.
#' @param observedPeak A list with names 'start' and 'end' containing
#' scalar values representing the calculated peak boundary points
#' @param filename A string containing the name of the current data file being
#' analyzed.
#'
#' @details A distribution is created from the set of all ppm values identified.
#' The most dense peak of this distribution is assumed to represent the standard
#' ppm error of the data.
#'
#' @return This function returns a scalar value representing ppm error estimate.
filterPpmError <- function(approvedPeaks, useGap, varExpThresh,
                            returnPpmPlots, plotDir, observedPeak,
                            filename) {

    ppmObs <- approvedPeaks$meanPPM
    ppmObs <- strsplit(split = ";", x = as.character(ppmObs))
    ppmObs <- lapply(ppmObs, as.numeric)
    ppmObs <- unlist(ppmObs)

    ## 2019-06-19
    ## corner case when all error measurements are identical.
    if(diff(range(ppmObs)) < .Machine$double.eps ^ 0.5) {
        stop(paste("All calculated ppm values are identical.",
                    "Error of data may be higher than the mass threshold value."
                    ))
    }

    message("-------- Number of ppm value across bins: ", length(ppmObs))
    if(length(ppmObs) > 10000) {
        ppmObs <- ppmObs[sample(x = seq_along(ppmObs), size = 5000)]
    }

    if(length(ppmObs) > 750) {

        checkPpm <- length(ppmObs)/2
        subsample <- TRUE
        while(subsample) {

            origDist <- stats::density(ppmObs, bw = 1)$y
            newDist1 <-  stats::density(sample(ppmObs, checkPpm), bw = 1)$y
            newDist2 <-  stats::density(sample(ppmObs, checkPpm), bw = 1)$y
            newDist3 <-  stats::density(sample(ppmObs, checkPpm), bw = 1)$y
            newDist4 <-  stats::density(sample(ppmObs, checkPpm), bw = 1)$y
            newDist5 <-  stats::density(sample(ppmObs, checkPpm), bw = 1)$y
            newDist6 <-  stats::density(sample(ppmObs, checkPpm), bw = 1)$y
            newDist7 <-  stats::density(sample(ppmObs, checkPpm), bw = 1)$y

            klDistance <- list()
            subSamples <- ls()[grep("newDist",ls())]
            for(j in seq_along(subSamples)) {
                klDistance[[j]] <- suppressWarnings(entropy::KL.empirical(
                    origDist,get(subSamples[j])))
            }
            klDistance <- unlist(klDistance)

            if(any(klDistance >= 0.5)) {
                subsample <- FALSE
            } else {
                checkPpm <- checkPpm/2
            }

        }

        ppmObs <- sample(ppmObs, checkPpm)
        message("-------- Number of ppm value across bins after",
                " KL Distance Filtering: ",
                length(ppmObs))

    }


    ## 2019-04-09 added this here since it doesn't make sense to cluster too few
    ## features
    if(length(ppmObs) < 100) {

        kmeansPPM <- kmeans(ppmObs, 1)

    } else if(useGap) {

        gapStat <- cluster::clusGap(x = as.matrix(ppmObs),
                                    FUNcluster = kmeans,
                                    K.max = 5,
                                    B = 7,
                                    verbose = FALSE)

        gapStat <- gapStat$Tab
        gap <- diff(-gapStat[,3]) > 0
        if(any(gap)) {
            clusters <- max(which(gap)) + 1
        } else {
            clusters <- 1
        }

        kmeansPPM <- kmeans(ppmObs, clusters)

    } else {


        ## estimating clustering based on hard coded 80% Vexp threshold
        clustCount <- 1
        varExp <- 0
        while(varExp < varExpThresh && clustCount < length(ppmObs)/2) {
            kmeansPPM <- kmeans(ppmObs, clustCount)
            varExp <- kmeansPPM$betweenss/kmeansPPM$totss
            clustCount <- clustCount + 1
        }

    }

    ## cluster which contains smallest ppm values
    clusterSize <- table(kmeansPPM$cluster)
    clusterSize <- sort(clusterSize, decreasing = TRUE)

    maxCluster <- names(clusterSize)[1]
    minCluster <- which(kmeansPPM$cluster == maxCluster)
    rm(clusterSize)

    x <- ppmObs
    n <- length(x)
    h <- 1
    ## delete this later
    gauss <- function(x) 1/sqrt(2*pi) * exp(-(x^2)/2)
    gaussDKE <- function(a, x) gauss((x - a)/h)/(n * h)

    bumps <- vapply(X = ppmObs[minCluster], FUN = gaussDKE, x = x,
                    FUN.VALUE = numeric(length = length(x)))
    wholeKDE <- vapply(X = ppmObs, FUN = gaussDKE,
                       x,
                       FUN.VALUE = numeric(length =
                                               length(x)))

    ## calculating this ahead of time to avoid unnecessary downstream
    ## math
    cKdeMean <- sum(rowSums(bumps))/length(minCluster)
    OutlierScore <- rowSums(wholeKDE)/(cKdeMean)

    scoreSub <- which(OutlierScore > 1)


    ppmEst <- max(ppmObs[scoreSub])
    maxX <- ppmEst
    ppmEst <- ppmEst + sd(ppmObs[scoreSub])*3


    if(returnPpmPlots) {


        title <- paste(filename,
                       'ppm Distribution of Bounded Peak\n Range (s):',
                        signif(observedPeak$start, digits = 4),
                        "-",
                        signif(observedPeak$end, digits = 4))
        title <- trimws(title)

        output <- file.path(plotDir,paste0(gsub(" ", "_", title), ".pdf"))
        output <- sub(":", "", output)

        ## error here...
        par(mar=c(1,1,1,1))
        grDevices::pdf(output, width = 8, height = 6)

        ## adding heuristic here to make ploting easier to see
        if(length(ppmObs) < 300) {
            bw <- 0.05
        } else {
            bw <- .1
        }

        plot(stats::density(ppmObs,bw = bw),
            main = title,
            cex.main = 1.2, cex.lab = 1.3, cex.axis = 1.2,
            xlab = "ppm Values") #+
        abline(v = maxX, lty = 2, col = "red") +
        abline(v = ppmEst, lty = 3, col = "blue")
        legend("topright",
                legend = c(paste("Maximum Outlier Score > 1:", signif(maxX,digits = 3)),
                    paste("ppm Estimate:", signif(ppmEst,digits = 3))),
                col = c("red","blue"),
                lty = c(2,3),cex = 1.1)
        grDevices::dev.off()

    }

    return(ppmEst)
}
