context("Preparing things for Autotuner")

library(mtbls2)


rawPaths <- c(
    system.file("mzData/MSpos-Ex2-cyp79-48h-Ag-1_1-B,3_01_9828.mzData",
                package = "mtbls2"),
    system.file("mzData/MSpos-Ex2-cyp79-48h-Ag-2_1-B,4_01_9830.mzData",
                package = "mtbls2"),
    system.file("mzData/MSpos-Ex2-cyp79-48h-Ag-4_1-B,4_01_9834.mzData",
                package = "mtbls2"))


runfile <- read.table(system.file(
    "a_mtbl2_metabolite_profiling_mass_spectrometry.txt", package = "mtbls2"),
    header = TRUE, stringsAsFactors = FALSE)

runfile <- runfile[sub("mzData/", "", runfile$Raw.Spectral.Data.File) %in%
                         basename(rawPaths),]

## Loading Autotuner
Autotuner <- Autotuner::createAutotuner(rawPaths,
                                        runfile,
                                        file_col = "Raw.Spectral.Data.File",
                                        factorCol = "Factor.Value.genotype.")

test_that(desc = "Auotuner object creaction",
          code = {

              # this function runs all tests when called upon

              ## Checking if Autotuner is able to correctly subset data
              ## by metadata columns
              expect_equal(sum(c(Autotuner@file_col,
                                 Autotuner@factorCol) %in%
                                    colnames(Autotuner@metadata)), 2)

              ## checking that Autotuner object is constructed with correct
              ## output

              expect_equal(class(Autotuner@time), "list")
              expect_equal(class(unlist(Autotuner@time)), "numeric")
              expect_equal(class(Autotuner@intensity), "list")
              expect_equal(class(unlist(Autotuner@intensity)), "numeric")
              expect_equal(class(Autotuner@file_paths), "character")
              expect_equal(class(unlist(Autotuner@intensity)), "numeric")

              #Autotuner@metadata
              #Autotuner@file_paths


})

lag <- 20
threshold<- 3
influence <- 0.1
signal <- lapply(getAutoIntensity(Autotuner), ThresholdingAlgo, lag, threshold,
                 influence)

test_that(desc = "Signal Processing Structure", code = {

    ## check that the correct object is returned from signal function
    expect_equal(length(signal), 3)
    expect_equal(length(signal[[1]]), 3)
    expect_equal(length(signal[[2]]), 3)
    expect_equal(length(signal[[3]]), 3)

})



test_that(desc = "Signal Processing Output",
          code = {

              ## check that computation took place
              naCheck <- list()
              for(i in seq_along(signal)) {
                  naCheck[[i]] <- sum(sapply(signal[[1]], function(x) {
                      all(is.na(x))
                  }))

              }

              expect_equal(naCheck[[1]], 0)
              expect_equal(naCheck[[2]], 0)
              expect_equal(naCheck[[3]], 0)

})

p <- plot_signals(Autotuner,
             threshold,
             ## index for which data files should be displayed
             sample_index = 1:3,
             signals = signal)

test_that(desc = "Checking plot signals function",
          code = {
              expect_equal(all(names(p$rect) == c("w", "h", "left", "top")),
                           TRUE)
          })


Autotuner <- isolatePeaks(Autotuner, returned_peaks = 10, signal)

test_that(desc = "Checking Function to Return Peaks",
          code = {

              nullCount <- sum(sapply(Autotuner@peaks, is.null))
              expect_equal(nullCount, 0)
              expect_equal(ncol(Autotuner@peaks[[1]]) <= 10, TRUE)
              expect_equal(length(unique(Autotuner@peak_difference$index)) > 3, TRUE)

})

p <- plot_peaks(Autotuner = Autotuner,
           boundary = 100,
           peak = 1)


test_that(desc = "Checking that plot peaks function works",
          code = {
              expect_equal(all(names(p$rect) == c("w", "h", "left", "top")),
                           TRUE)
          })


test_that(desc = "Checking Peakwidth_table",
          code = {
                expect_equal(class(Autotuner@peak_table), "data.frame")
                expect_equal(any(is.na(Autotuner@peak_table)), FALSE)
          })

test_that(desc = "Checking peak_time_difference",
          code = {
              expect_equal(class(Autotuner@peak_difference), "data.frame")
          })


eicParamEsts <- EICparams(Autotuner = Autotuner,
                          massThresh = .005,
                          verbose = FALSE,
                          returnPpmPlots = FALSE,
                          useGap = TRUE)


test_that(desc = "checking EICparams function",
          code = {
              expect_equal(class(eicParamEsts), "data.frame")
              expect_equal(length(unique(eicParamEsts$sampleID)), 3)
              expect_equal(max(eicParamEsts$ppm) > 5, TRUE)
              expect_equal(max(eicParamEsts$maxPw) < 30, TRUE)
          })
