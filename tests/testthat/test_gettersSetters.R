context("Testing Accessor Functions")

Autotuner <- readRDS(system.file("extdata/AutotunerObj.rds", package="Autotuner"))

intensity <- getAutoIntensity(Autotuner)
time <- getAutoTime(Autotuner)
peaks <- getAutoPeaks(Autotuner)
peak_table <- getAutoPeak_table(Autotuner)
peak_difference <- getAutoPeak_difference(Autotuner)
metadata <- getAutoMetadata(Autotuner)
file_paths <- getAutoFile_paths(Autotuner)
file_col <- getAutoFile_col(Autotuner)
factorCol <- getAutoFactorCol(Autotuner)

test_that(desc = "Auotuner getters check",
          code = {
              expect_equal(length(intensity), 3)
              expect_equal(all(sapply(intensity, class) == "numeric"), TRUE)
              expect_equal(length(time), 3)
              expect_equal(all(sapply(time, class) == "numeric"), TRUE)
              expect_equal(length(peaks), 3)
              expect_equal(all(sapply(peaks, class) == "data.frame"), TRUE)
              expect_equal(nrow(peak_difference), 15)
              expect_equal(is(peak_difference)[1], "data.frame")
              expect_equal(nrow(metadata), 3)
              expect_equal(is(metadata)[1], "data.frame")
              expect_equal(length(file_paths), 3)
              expect_equal(is(file_paths)[1], "character")
              expect_equal(length(file_col), 1)
              expect_equal(is(file_col)[1], "character")
              expect_equal(length(factorCol), 1)
              expect_equal(is(factorCol)[1], "character")
          })


Autotuner <- setAutoIntensity(intensity = list(), Autotuner)
Autotuner <- setAutoTime(time = list(), Autotuner)
Autotuner <- setAutoPeaks(peaks = list(), Autotuner)
Autotuner <- setAutoPeak_table(peak_table = data.frame(), Autotuner)
Autotuner <- setAutoPeak_difference(peak_difference = data.frame(), Autotuner)
Autotuner <- setAutoMetadata(metadata = data.frame(), Autotuner)
Autotuner <- setAutoFile_paths(file_paths = character(), Autotuner)
Autotuner <- setAutoFile_col(file_col = character(), Autotuner)
Autotuner <- setAutoFactorCol(factorCol = character(), Autotuner)

test_that(desc = "Autotuner setter check",
          code = {
              expect_equal(length(getAutoIntensity(Autotuner)), 0)
              expect_equal(length(getAutoFactorCol(Autotuner)), 0)
              expect_equal(length(getAutoFile_col(Autotuner)), 0)
              expect_equal(length(getAutoFile_paths(Autotuner)), 0)
              expect_equal(length(getAutoMetadata(Autotuner)), 0)
              expect_equal(length(getAutoPeak_difference(Autotuner)), 0)
              expect_equal(length(getAutoPeak_table(Autotuner)), 0)
              expect_equal(length(getAutoPeaks(Autotuner)), 0)
              expect_equal(length(getAutoTime(Autotuner)), 0)

})


