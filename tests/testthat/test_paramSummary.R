context("Performing dataset wide parameter selection")

Autotuner <- readRDS(system.file("data/Autotuner.rds", package="Autotuner"))
eicTable <- readRDS(system.file('data/eicParamEsts.rds', package = "Autotuner"))

returnParams <- returnParams(eicTable, Autotuner)

test_that(desc = "Testing Max Peakwidth Estimate",
          code = {
              expect_equal(class(returnParams), "list")
              expect_equal(length(returnParams), 2)
              expect_equal(dim(returnParams$ticParams), c(3,2))
              expect_equal(dim(returnParams$eicParams), c(7,4))
          })
