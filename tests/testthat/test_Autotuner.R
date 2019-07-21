context("Performing Parameter Selection by Autotuner")


## loading in the files to test if AutoTuner is correct
useGap <- TRUE
returnPpmPlots <- FALSE
massThresh <- 0.005
filename <- ''
data("mzDb", package="Autotuner")
header <- readRDS(system.file("extdata/header.rds",
                              package = "Autotuner"))
data("observedPeak", package="Autotuner")


scanDiff <- diff(header$retentionTime[header$msLevel == 1L])
rate <- mean(scanDiff)
rm(scanDiff)

sortedAllEIC <- dissectScans(mzDb,
                             observedPeak = observedPeak,
                             header = header)

test_that(desc = "Checking dissectScans",
          code = {
              expect_false(nrow(sortedAllEIC) == 0)
              expect_equal(sum(is.na(sortedAllEIC)), 0)
              expect_equal(unique(apply(sortedAllEIC, 2, class)), "numeric")
          })

boundaries <- range(sortedAllEIC$scanID)

# Checking if sorted scans pass mz error threshold ------------------------
matchedMasses <- rle(diff(sortedAllEIC$mz) < massThresh)

noiseAndPeaks <- filterPeaksfromNoise(matchedMasses)
no_match <- noiseAndPeaks[[1]]
truePeaks <- noiseAndPeaks[[2]]
rm(noiseAndPeaks)

test_that(desc = "Testing noise filtering",
          code = {
              expect_false(length(no_match) == 0)
              expect_false(length(truePeaks) == 0)
          })

approvedPeaks <- findTruePeaks(truePeaks, sortedAllEIC)

test_that(desc = "Testing findTruePeaks",
          code = {
              expect_equal(class(approvedPeaks), "data.frame")
              expect_false(nrow(approvedPeaks) == 0)
          })

overlappingScans <- sum(approvedPeaks$multipleInScan)

ppmEst <- suppressMessages(filterPpmError(approvedPeaks, useGap, varExpThresh,
                         returnPpmPlots, plotDir, observedPeak,
                         filename))

test_that(desc = "Testing ppm estimation",
          code = {
              expect_equal(length(ppmEst), 1)
              expect_equal(class(ppmEst), "numeric")
              expect_false(is.na(ppmEst))
          })

ppmObs <- approvedPeaks$meanPPM
ppmObs <- strsplit(split = ";", x = as.character(ppmObs))
ppmObs <- sapply(ppmObs, as.numeric)


noisyBin <- lapply(ppmObs, function(ppm) {
    any(ppm > ppmEst)
})
noisyBin <- unlist(noisyBin)
approvScorePeaks <- approvedPeaks[!noisyBin,]
SNest <- estimateSNThresh(no_match,
                          sortedAllEIC, approvScorePeaks)

test_that(desc = "Testing s/n Estimation",
          code = {
              expect_equal(class(SNest), "numeric")
              expect_false(any(is.na(SNest)))
          })

SNest <- min(SNest)

scanEst <- min(approvScorePeaks$scanCount)

### Noise Intensity Estimate
noiseEst <- min(approvScorePeaks$minIntensity) - 100
if(noiseEst < 0) {
    noiseEst <- min(approvScorePeaks$minIntensity) + 10
}

### Prefilter Intensity Estimate
intensityEst <- min(approvScorePeaks$Intensity)/sqrt(2)

### peakWidth Estimate
boundaries <- range(sortedAllEIC$scanID)
maxPw <- findPeakWidth(approvScorePeaks = approvScorePeaks,
                       mzDb = mzDb,
                       header = header,
                       sortedAllEIC = sortedAllEIC,
                       boundaries = boundaries,
                       ppmEst = ppmEst)

test_that(desc = "Testing Max Peakwidth Estimate",
          code = {
              expect_equal(class(maxPw), "numeric")
              expect_false(is.na(maxPw))
          })

